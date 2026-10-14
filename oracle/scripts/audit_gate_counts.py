"""Reconcile component counts with the full-oracle CCX totals."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
paths = {d: ROOT / "results/gate-audit" / f"{d}.json" for d in ["bonnetain", "revised"]}
records = {d: json.loads(f.read_text()) for d, f in paths.items()}
for d, a in records.items():
    assert a["gate_audit"]["noninteger_append_weights"] == []
    for name, layer in a["gate_audit"]["layers"].items():
        assert not layer["unresolved_direct_names"], (d, name)
        assert layer["known_direct_ccx_sum"] == layer["reported_ccx"], (d, name)
    assert a["resources"]["ccx"] < 2**53
    published = ROOT / (
        "results/bonnetain/oracle.json"
        if d == "bonnetain"
        else "results/reference/oracle.json"
    )
    assert a["resources"] == json.loads(published.read_text())["resources"]

B, R = records["bonnetain"], records["revised"]


def cc(a, name):
    return a["layer_resources"][name]["ccx"]


Cb = cc(B, "BonnetainPositiveCharacterFlag")
Cr = cc(R, "UnitaryPositiveCharacterFlag")
char_saving = 4 * (Cb - Cr)
local_increase = cc(R, "UnitaryValueImageFromPoint") - cc(
    B, "UnitaryValueImageFromPoint"
)
assert cc(B, "UnitaryArithmeticPrepare") - cc(R, "UnitaryArithmeticPrepare") == 2 * (
    Cb - Cr
)
assert cc(B, "UnitaryArithmeticUnprepare") - cc(
    R, "UnitaryArithmeticUnprepare"
) == 2 * (Cb - Cr)
for a in [B, R]:
    assert cc(a, "UnitaryValueOOP") == sum(
        cc(a, n)
        for n in [
            "UnitaryArithmeticPrepare",
            "UnitaryValueImageFromPoint",
            "UnitaryArithmeticUnprepare",
        ]
    )
    assert cc(a, "ValueSignedInplace") == 2 * cc(a, "ValueSignedOOP")
    assert a["resources"]["ccx"] == 1300 * (
        cc(a, "ValueSignedInplace") + 2 * cc(a, "UnitaryActionSelector")
    )
net = local_increase - char_saving
assert net == cc(R, "UnitaryValueOOP") - cc(B, "UnitaryValueOOP")
assert 2600 * net == R["resources"]["ccx"] - B["resources"]["ccx"]
assert Cb == 1924 * 8625671 + 2039
assert Cr == 2300 * 2874883 + 2039
point_rows = []
for a, name in [(B, "BonnetainPointStep"), (R, "ExplicitUnitaryWPointStep")]:
    layer = a["gate_audit"]["layers"][name]
    operations = [
        e
        for e in layer["direct_children"]
        if any(s in e["class"] for s in ["ModularProductAdd", "ModularSquareAdd"])
    ]
    subtotal = sum(int(e["quantity"]) * e["ccx_each"] for e in operations)
    point_rows.append(
        dict(
            design=a["design"],
            component=name,
            operations=operations,
            arithmetic_calls=sum(int(e["quantity"]) for e in operations),
            multiplicative_ccx=subtotal,
            other_ccx=layer["reported_ccx"] - subtotal,
            total_ccx=layer["reported_ccx"],
        )
    )
result = dict(
    outer_calls=dict(
        M=1300, signed_oop_per_inplace=2, total_oop=2600, character_calls_per_oop=4
    ),
    character=dict(
        bonnetain_ccx=Cb,
        revised_ccx=Cr,
        bonnetain_bit_steps=1924,
        bonnetain_field_calls=5772,
        revised_field_calls=2300,
        shared_equality_overhead_ccx=2039,
        native_reduction_percent=100 * (1 - Cr / Cb),
    ),
    per_oop=dict(
        character_ccx_saved=char_saving,
        local_mapping_ccx_increase=local_increase,
        net_ccx_increase=net,
    ),
    full=dict(ccx_increase=2600 * net, T7_increase=7 * 2600 * net),
    point_steps=point_rows,
)
output = ROOT / "results/gate-audit/summary.json"
output.write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        {
            k: result[k]
            for k in [
                "outer_calls",
                "character",
                "per_oop",
                "full",
            ]
        },
        indent=2,
    )
)
