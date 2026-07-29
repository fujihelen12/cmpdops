import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import math
import string
import re
from io import BytesIO


# ==================================================
# PAGE CONFIG
# ==================================================

st.set_page_config(
    page_title="CompoundOps Assistant",
    layout="wide"
)

st.title("CompoundOps Assistant")

with st.sidebar:

    st.header("Application Controls")

    if st.button(
        "🗑️ Clear Entire Session",
        use_container_width=True
    ):
        for key in list(st.session_state.keys()):
            del st.session_state[key]

        st.rerun()
st.write(
    "Inventory intelligence, assay design, solvent-normalized backfill planning, "
    "intermediate dilution strategy, and dose-response plate mapping."
)


# ==================================================
# LOAD INVENTORY
# ==================================================

inventory = pd.read_csv(
    "inventory.csv",
    encoding="cp1252"
)


# ==================================================
# STANDARDIZE COLUMN NAMES
# ==================================================

inventory = inventory.rename(
    columns={
        "Compound ID": "compound_id",
        "MW": "mw",
        "Stock (mM)": "stock_mM",
        "Available (µL)": "available_uL",
        "Freeze-Thaw": "freeze_thaw",
        "Status": "inventory_status",
        "Solvent": "solvent",
        "LogP": "logP",
        "Storage Temp": "storage_temp",
        "Container Type": "container_type",
        "Purity %": "purity_percent",
        "Expiry": "expiry",
        "Location": "location"
    }
)


# ==================================================
# REQUIRED COLUMN SAFETY
# ==================================================

required_cols = [
    "compound_id",
    "mw",
    "stock_mM",
    "available_uL",
    "freeze_thaw",
    "inventory_status",
    "solvent",
    "logP",
    "storage_temp",
    "container_type",
    "purity_percent",
    "expiry",
    "location"
]

for col in required_cols:
    if col not in inventory.columns:
        inventory[col] = ""


# ==================================================
# DATA CLEANUP
# ==================================================

inventory["compound_id"] = (
    inventory["compound_id"]
    .astype(str)
    .str.strip()
)

inventory["location"] = (
    inventory["location"]
    .astype(str)
    .str.strip()
)

inventory["container_type"] = (
    inventory["container_type"]
    .astype(str)
    .str.strip()
)

inventory["solvent"] = (
    inventory["solvent"]
    .astype(str)
    .str.strip()
)

inventory["inventory_status"] = (
    inventory["inventory_status"]
    .astype(str)
    .str.strip()
)

numeric_cols = [
    "mw",
    "stock_mM",
    "available_uL",
    "freeze_thaw",
    "logP",
    "purity_percent"
]

for col in numeric_cols:
    inventory[col] = pd.to_numeric(
        inventory[col],
        errors="coerce"
    )

inventory = inventory.reset_index(drop=True)

inventory["Aliquot ID"] = (
    inventory["compound_id"].astype(str)
    + " | "
    + inventory["location"].astype(str)
    + " | "
    + inventory["container_type"].astype(str)
    + " | Row "
    + (inventory.index + 1).astype(str)
)


# ==================================================
# GENERAL HELPERS
# ==================================================

def clean_unique_values(series):
    values = (
        series
        .dropna()
        .astype(str)
        .str.strip()
    )

    values = [
        v for v in values.unique().tolist()
        if v and v.lower() != "nan"
    ]

    return sorted(values)


def parse_compound_ids(text):
    if not text:
        return []

    cleaned = (
        text.replace(",", "\n")
        .replace(";", "\n")
        .replace("\t", "\n")
        .replace(" ", "\n")
        .splitlines()
    )

    ids = []

    for item in cleaned:
        item = item.strip()

        if item:
            ids.append(item)

    return ids


def convert_concentration_to_mM(value, unit):
    value = float(value)

    if unit == "mM":
        return value

    if unit == "µM":
        return value / 1000

    if unit == "nM":
        return value / 1_000_000

    if unit == "pM":
        return value / 1_000_000_000

    return value


def format_concentration(value_mM):
    if value_mM is None or pd.isna(value_mM):
        return ""

    value_mM = float(value_mM)

    if value_mM >= 1:
        return f"{value_mM:.3g} mM"

    value_uM = value_mM * 1000

    if value_uM >= 1:
        return f"{value_uM:.3g} µM"

    value_nM = value_mM * 1_000_000

    if value_nM >= 1:
        return f"{value_nM:.3g} nM"

    value_pM = value_mM * 1_000_000_000

    return f"{value_pM:.3g} pM"


def format_volume_nL(value_nL):
    if value_nL is None or pd.isna(value_nL):
        return ""

    value_nL = float(value_nL)

    if value_nL >= 1:
        return f"{value_nL:.3g} nL"

    value_pL = value_nL * 1000

    return f"{value_pL:.3g} pL"


def format_volume_uL(value_uL):
    if value_uL is None or pd.isna(value_uL):
        return ""

    return f"{float(value_uL):.3g} µL"


def build_concentration_series(
    top_conc_mM,
    dose_points,
    dilution_factor
):
    concentrations = []
    current = float(top_conc_mM)

    for _ in range(int(dose_points)):
        concentrations.append(current)
        current = current / float(dilution_factor)

    return concentrations


def calculate_solvent_backfill(
    transfer_nL,
    final_assay_volume_uL,
    target_solvent_percent,
    source_solvent_percent=100
):
    if transfer_nL is None or pd.isna(transfer_nL):
        return {
            "Target Total Solvent (nL)": None,
            "Solvent From Transfer (nL)": None,
            "Solvent Backfill (nL)": None,
            "Aqueous Backfill (nL)": None,
            "Actual Final Solvent %": None,
            "Solvent Backfill Status": "Review"
        }

    final_assay_volume_nL = final_assay_volume_uL * 1000

    target_total_solvent_nL = (
        final_assay_volume_nL
        * target_solvent_percent
        / 100
    )

    solvent_from_transfer_nL = (
        transfer_nL
        * source_solvent_percent
        / 100
    )

    solvent_backfill_nL = (
        target_total_solvent_nL
        - solvent_from_transfer_nL
    )

    status = "OK"

    if solvent_backfill_nL < 0:
        solvent_backfill_nL = 0
        status = "Transfer exceeds target solvent"

    aqueous_backfill_nL = (
        final_assay_volume_nL
        - transfer_nL
        - solvent_backfill_nL
    )

    if aqueous_backfill_nL < 0:
        aqueous_backfill_nL = 0
        status = "Total volume exceeded"

    actual_final_solvent_percent = (
        (
            solvent_from_transfer_nL
            + solvent_backfill_nL
        )
        / final_assay_volume_nL
        * 100
    )

    return {
        "Target Total Solvent (nL)": target_total_solvent_nL,
        "Solvent From Transfer (nL)": solvent_from_transfer_nL,
        "Solvent Backfill (nL)": solvent_backfill_nL,
        "Aqueous Backfill (nL)": aqueous_backfill_nL,
        "Actual Final Solvent %": actual_final_solvent_percent,
        "Solvent Backfill Status": status
    }


# ==================================================
# INVENTORY FUNCTIONS
# ==================================================

def inventory_health(row):
    flags = []

    status = str(row["inventory_status"]).strip()

    if status != "Available":
        flags.append("Inventory unavailable")

    if pd.notna(row["freeze_thaw"]):
        if row["freeze_thaw"] > 4:
            flags.append("Freeze-thaw exceeds SOP guideline")

    expiry = pd.to_datetime(
        row["expiry"],
        errors="coerce"
    )

    if pd.notna(expiry):
        if expiry < pd.Timestamp.today().normalize():
            flags.append("Expired")

    if pd.notna(row["available_uL"]):
        usable = row["available_uL"] - 10

        if usable <= 0:
            flags.append("Low usable volume")

    location = str(row["location"]).strip()

    if location == "" or location.lower() == "nan":
        flags.append("Location missing")

    if len(flags) == 0:
        return "Healthy"

    return " | ".join(flags)


def solvent_guidance(solvent):
    solvent = str(solvent).upper()

    if solvent == "DMSO":
        return """
**Storage**
- Store frozen when possible.

**Handling**
- Minimize moisture exposure.
- Inspect for precipitation before use.

**Special Considerations**
- Final assay DMSO should remain within assay limit.
- Vehicle controls are recommended.
"""

    if solvent == "METHANOL":
        return """
**Storage**
- Keep tightly sealed.

**Handling**
- Minimize evaporation.
- Inspect before use.

**Special Considerations**
- Review assay compatibility before use.
"""

    if solvent == "ETHANOL":
        return """
**Storage**
- Keep tightly sealed.

**Handling**
- Minimize evaporation.

**Special Considerations**
- Review assay tolerance before use.
"""

    if solvent == "ACETONITRILE":
        return """
**Storage**
- Follow documented storage conditions.

**Handling**
- Minimize open exposure.

**Special Considerations**
- Review assay compatibility before use.
"""

    if solvent in [
        "PBS",
        "BUFFER",
        "WATER"
    ]:
        return """
**Storage**
- Verify storage temperature.

**Handling**
- Inspect before use.

**Special Considerations**
- Review stability before use.
"""

    return """
**SOP Guidance**
- Review compound-specific handling requirements.
"""


def run_inventory_global_search(
    inventory_df,
    search_term,
    solvent_filter,
    container_filter,
    status_filter,
    location_filter
):
    filtered = inventory_df.copy()

    search_term = str(search_term).strip()

    if search_term:
        mask = (
            filtered
            .astype(str)
            .apply(
                lambda col: col.str.contains(
                    search_term,
                    case=False,
                    na=False
                )
            )
            .any(axis=1)
        )

        filtered = filtered[mask].copy()

    if solvent_filter:
        filtered = filtered[
            filtered["solvent"].isin(solvent_filter)
        ].copy()

    if container_filter:
        filtered = filtered[
            filtered["container_type"].isin(container_filter)
        ].copy()

    if status_filter:
        filtered = filtered[
            filtered["inventory_status"].isin(status_filter)
        ].copy()

    if location_filter:
        filtered = filtered[
            filtered["location"].isin(location_filter)
        ].copy()

    return filtered


# ==================================================
# ASSAY DESIGN FUNCTIONS
# ==================================================

def direct_assay_row_recommendation(reasons):
    if len(reasons) == 0:
        return "Direct stock transfer appears suitable for this dose point."

    if any("below Echo minimum" in r for r in reasons):
        return "This point likely requires an intermediate dilution because the direct transfer is below the Echo minimum."

    if any("above Echo maximum" in r for r in reasons):
        return "This point may require a stronger source or revised assay volume because transfer exceeds the Echo maximum."

    if any("exceeds stock concentration" in r for r in reasons):
        return "Requested target concentration exceeds source stock concentration. Lower the top dose or locate a stronger stock."

    if any("target solvent" in r for r in reasons):
        return "Compound transfer alone exceeds the target solvent percentage. Lower dose, use stronger stock, or adjust assay conditions."

    if any("inventory" in r.lower() for r in reasons):
        return "Review inventory record before using this compound in assay design."

    return "Review assay design parameters before proceeding."


def build_direct_stock_assay_design(
    selected_compounds,
    concentrations_mM,
    assay_well_volume_uL,
    min_transfer_nL,
    max_transfer_nL,
    plates,
    replicates,
    target_solvent_percent,
    source_solvent_percent
):
    rows = []
    final_assay_volume_nL = assay_well_volume_uL * 1000

    for _, compound in selected_compounds.iterrows():

        compound_id = compound["compound_id"]
        stock_mM = compound["stock_mM"]
        available_uL = compound["available_uL"]
        usable_uL = available_uL - 10 if pd.notna(available_uL) else None
        solvent = str(compound["solvent"])
        status = str(compound["inventory_status"])
        freeze_thaw = compound["freeze_thaw"]
        location = compound["location"]
        aliquot_id = compound.get("Aliquot ID", f"{compound_id} | {location}")

        expiry = pd.to_datetime(
            compound["expiry"],
            errors="coerce"
        )

        for idx, target_mM in enumerate(
            concentrations_mM,
            start=1
        ):

            reasons = []

            if pd.isna(stock_mM) or stock_mM <= 0:
                transfer_nL = None
                reasons.append("Missing or invalid stock concentration")

            else:
                transfer_nL = (
                    target_mM
                    / stock_mM
                    * final_assay_volume_nL
                )

            if pd.notna(stock_mM):
                if target_mM > stock_mM:
                    reasons.append(
                        "Target concentration exceeds stock concentration"
                    )

            if transfer_nL is not None:

                if transfer_nL < min_transfer_nL:
                    reasons.append(
                        "Direct transfer below Echo minimum"
                    )

                if transfer_nL > max_transfer_nL:
                    reasons.append(
                        "Direct transfer above Echo maximum"
                    )

                required_source_uL = (
                    transfer_nL
                    * int(plates)
                    * int(replicates)
                    / 1000
                )

            else:
                required_source_uL = None

            if usable_uL is not None and required_source_uL is not None:
                if required_source_uL > usable_uL:
                    reasons.append(
                        "Required source volume exceeds usable inventory"
                    )

            if status != "Available":
                reasons.append(
                    "Inventory status is not Available"
                )

            if pd.notna(freeze_thaw):
                if freeze_thaw > 4:
                    reasons.append(
                        "Freeze-thaw exceeds SOP guideline"
                    )

            if pd.notna(expiry):
                if expiry < pd.Timestamp.today().normalize():
                    reasons.append(
                        "Compound is expired"
                    )

            backfill = calculate_solvent_backfill(
                transfer_nL=transfer_nL,
                final_assay_volume_uL=assay_well_volume_uL,
                target_solvent_percent=target_solvent_percent,
                source_solvent_percent=source_solvent_percent
            )

            if backfill["Solvent Backfill Status"] != "OK":
                reasons.append(
                    "Compound transfer exceeds target solvent normalization"
                )

            if len(reasons) == 0:
                assay_status = "Direct OK"
            else:
                assay_status = "Review"

            rows.append(
                {
                    "Compound ID": compound_id,
                    "Aliquot ID": aliquot_id,
                    "Dose Point": idx,
                    "Target Conc (mM)": target_mM,
                    "Target Conc": format_concentration(target_mM),
                    "Stock (mM)": stock_mM,
                    "Solvent": solvent,
                    "Source Solvent %": source_solvent_percent,
                    "Assay Well Volume (uL)": assay_well_volume_uL,
                    "Transfer Numeric (nL)": transfer_nL,
                    "Compound Transfer": format_volume_nL(transfer_nL),
                    "Target Total Solvent": format_volume_nL(
                        backfill["Target Total Solvent (nL)"]
                    ),
                    "Solvent From Compound": format_volume_nL(
                        backfill["Solvent From Transfer (nL)"]
                    ),
                    "Solvent Backfill": format_volume_nL(
                        backfill["Solvent Backfill (nL)"]
                    ),
                    "Aqueous / Media Backfill": format_volume_nL(
                        backfill["Aqueous Backfill (nL)"]
                    ),
                    "Actual Final Solvent %": backfill[
                        "Actual Final Solvent %"
                    ],
                    "Required Source Volume (uL)": required_source_uL,
                    "Required Source Volume": format_volume_uL(
                        required_source_uL
                    ),
                    "Usable Inventory (uL)": usable_uL,
                    "Assay Status": assay_status,
                    "Review Reasons": (
                        " | ".join(reasons)
                        if reasons
                        else "None"
                    ),
                    "Recommendation": direct_assay_row_recommendation(
                        reasons
                    )
                }
            )

    return pd.DataFrame(rows)


# ==================================================
# INTERMEDIATE STRATEGY FUNCTIONS
# ==================================================

def group_indices_evenly(indices, group_count):
    groups = []

    if len(indices) == 0:
        return groups

    group_size = math.ceil(
        len(indices) / group_count
    )

    for group_number in range(group_count):

        group = indices[
            group_number * group_size:
            (group_number + 1) * group_size
        ]

        if group:
            groups.append(group)

    return groups


def calculate_intermediate_recipe(
    compound_id,
    stock_mM,
    solvent,
    group_number,
    group_indices,
    concentrations_mM,
    assay_well_volume_uL,
    preferred_intermediate_transfer_nL,
    intermediate_working_volume_uL,
    intermediate_dead_volume_uL,
    intermediate_solvent_percent,
    dilution_buffer_name
):
    final_assay_volume_nL = assay_well_volume_uL * 1000

    first_idx = group_indices[0]
    first_target_mM = concentrations_mM[first_idx]

    intermediate_conc_mM = (
        first_target_mM
        * final_assay_volume_nL
        / preferred_intermediate_transfer_nL
    )

    if pd.notna(stock_mM) and stock_mM > 0:
        intermediate_conc_mM = min(
            intermediate_conc_mM,
            stock_mM
        )

    intermediate_conc_mM = max(
        intermediate_conc_mM,
        first_target_mM
    )

    intermediate_id = f"{compound_id}_INT_{group_number}"

    total_intermediate_volume_uL = (
        intermediate_working_volume_uL
        + intermediate_dead_volume_uL
    )

    recipe_flags = []

    if pd.isna(stock_mM) or stock_mM <= 0:

        source_stock_uL = None
        solvent_backfill_uL = None
        diluent_uL = None
        actual_intermediate_solvent_percent = None
        recipe_flags.append(
            "Missing or invalid source stock concentration"
        )

    else:

        source_stock_uL = (
            intermediate_conc_mM
            / stock_mM
            * total_intermediate_volume_uL
        )

        target_solvent_volume_uL = (
            total_intermediate_volume_uL
            * intermediate_solvent_percent
            / 100
        )

        solvent_from_source_uL = source_stock_uL

        solvent_backfill_uL = (
            target_solvent_volume_uL
            - solvent_from_source_uL
        )

        if solvent_backfill_uL < 0:
            solvent_backfill_uL = 0
            recipe_flags.append(
                "Source stock volume already exceeds target intermediate solvent percent"
            )

        diluent_uL = (
            total_intermediate_volume_uL
            - source_stock_uL
            - solvent_backfill_uL
        )

        if diluent_uL < 0:
            diluent_uL = 0
            recipe_flags.append(
                "Intermediate recipe volume exceeds total intermediate volume"
            )

        actual_intermediate_solvent_percent = (
            (
                solvent_from_source_uL
                + solvent_backfill_uL
            )
            / total_intermediate_volume_uL
            * 100
        )

    dose_points_covered = [
        f"D{i + 1}"
        for i in group_indices
    ]

    return {
        "Compound ID": compound_id,
        "Intermediate ID": intermediate_id,
        "Intermediate Conc (mM)": intermediate_conc_mM,
        "Intermediate Conc": format_concentration(
            intermediate_conc_mM
        ),
        "Source Stock (mM)": stock_mM,
        "Solvent": solvent,
        "Intermediate Target Solvent %": intermediate_solvent_percent,
        "Actual Intermediate Solvent %": actual_intermediate_solvent_percent,
        "Working Volume (uL)": intermediate_working_volume_uL,
        "Dead Volume (uL)": intermediate_dead_volume_uL,
        "Total Intermediate Volume (uL)": total_intermediate_volume_uL,
        "Source Stock to Make Intermediate (uL)": source_stock_uL,
        "Matching Solvent Backfill to Make Intermediate (uL)": solvent_backfill_uL,
        "Diluent / Buffer": dilution_buffer_name,
        "Diluent Volume to Make Intermediate (uL)": diluent_uL,
        "Dose Points Covered": ", ".join(dose_points_covered),
        "Recipe Flags": (
            " | ".join(recipe_flags)
            if recipe_flags
            else "None"
        ),
        "Recipe Recommendation": (
            f"Prepare {intermediate_id} at "
            f"{format_concentration(intermediate_conc_mM)} "
            f"with approximately {intermediate_solvent_percent:.3g}% {solvent}. "
            f"This intermediate supports {', '.join(dose_points_covered)}."
        )
    }


def build_intermediate_strategy(
    selected_compounds,
    concentrations_mM,
    assay_well_volume_uL,
    min_transfer_nL,
    max_transfer_nL,
    target_final_solvent_percent,
    source_stock_solvent_percent,
    preferred_intermediate_transfer_nL,
    max_intermediates_per_compound,
    intermediate_working_volume_uL,
    intermediate_dead_volume_uL,
    intermediate_solvent_percent,
    dilution_buffer_name,
    intermediate_start_mode,
    manual_start_dose
):
    intermediate_rows = []
    assignment_rows = []

    final_assay_volume_nL = assay_well_volume_uL * 1000

    for _, compound in selected_compounds.iterrows():

        compound_id = compound["compound_id"]
        stock_mM = compound["stock_mM"]
        solvent = str(compound["solvent"])
        location = compound["location"]
        aliquot_id = compound.get("Aliquot ID", f"{compound_id} | {location}")

        direct_flags = []
        direct_transfer_list = []

        for target_mM in concentrations_mM:

            if pd.isna(stock_mM) or stock_mM <= 0:

                direct_transfer_nL = None
                can_direct = False

            else:

                direct_transfer_nL = (
                    target_mM
                    / stock_mM
                    * final_assay_volume_nL
                )

                can_direct = (
                    target_mM <= stock_mM
                    and direct_transfer_nL >= min_transfer_nL
                    and direct_transfer_nL <= max_transfer_nL
                )

            direct_transfer_list.append(
                direct_transfer_nL
            )

            direct_flags.append(
                can_direct
            )

        if intermediate_start_mode == "Auto from first non-direct dose":

            needs_intermediate = [
                idx
                for idx, flag in enumerate(direct_flags)
                if not flag
            ]

        else:

            start_idx = int(manual_start_dose) - 1

            needs_intermediate = [
                idx
                for idx in range(
                    start_idx,
                    len(concentrations_mM)
                )
            ]

        intermediate_groups = []

        if len(needs_intermediate) > 0:

            group_count = min(
                int(max_intermediates_per_compound),
                len(needs_intermediate)
            )

            intermediate_groups = group_indices_evenly(
                needs_intermediate,
                group_count
            )

        intermediate_lookup = {}

        for group_number, group_indices in enumerate(
            intermediate_groups,
            start=1
        ):

            recipe = calculate_intermediate_recipe(
                compound_id=compound_id,
                stock_mM=stock_mM,
                solvent=solvent,
                group_number=group_number,
                group_indices=group_indices,
                concentrations_mM=concentrations_mM,
                assay_well_volume_uL=assay_well_volume_uL,
                preferred_intermediate_transfer_nL=preferred_intermediate_transfer_nL,
                intermediate_working_volume_uL=intermediate_working_volume_uL,
                intermediate_dead_volume_uL=intermediate_dead_volume_uL,
                intermediate_solvent_percent=intermediate_solvent_percent,
                dilution_buffer_name=dilution_buffer_name
            )

            intermediate_rows.append(recipe)

            for idx in group_indices:

                intermediate_lookup[idx] = {
                    "Intermediate ID": recipe["Intermediate ID"],
                    "Intermediate Conc (mM)": recipe["Intermediate Conc (mM)"],
                    "Intermediate Solvent %": recipe["Actual Intermediate Solvent %"]
                }

        for idx, target_mM in enumerate(
            concentrations_mM,
            start=1
        ):

            zero_idx = idx - 1
            direct_transfer_nL = direct_transfer_list[zero_idx]

            source_type = "Review"
            source_id = "Unassigned"
            source_conc_mM = None
            transfer_nL = direct_transfer_nL
            source_solvent_percent_for_assay = source_stock_solvent_percent
            assignment_status = "Review"
            recommendation = "Review dose point."

            if (
                zero_idx not in intermediate_lookup
                and direct_flags[zero_idx]
            ):

                source_type = "Source Stock"
                source_id = "SOURCE"
                source_conc_mM = stock_mM
                transfer_nL = direct_transfer_nL
                source_solvent_percent_for_assay = source_stock_solvent_percent
                assignment_status = "Direct from Stock"
                recommendation = "Transfer directly from source stock."

            elif zero_idx in intermediate_lookup:

                inter = intermediate_lookup[zero_idx]

                source_type = "Intermediate"
                source_id = inter["Intermediate ID"]
                source_conc_mM = inter["Intermediate Conc (mM)"]

                if inter["Intermediate Solvent %"] is None:
                    source_solvent_percent_for_assay = intermediate_solvent_percent
                else:
                    source_solvent_percent_for_assay = inter["Intermediate Solvent %"]

                transfer_nL = (
                    target_mM
                    / source_conc_mM
                    * final_assay_volume_nL
                )

                if transfer_nL < min_transfer_nL:

                    assignment_status = "Low Transfer After Intermediate"

                    recommendation = (
                        "Intermediate assigned, but transfer is still below Echo minimum. "
                        "Consider increasing allowed intermediates or lowering intermediate concentration."
                    )

                elif transfer_nL > max_transfer_nL:

                    assignment_status = "High Transfer After Intermediate"

                    recommendation = (
                        "Intermediate assigned, but transfer exceeds Echo maximum. "
                        "Consider a stronger intermediate."
                    )

                else:

                    assignment_status = "From Intermediate"

                    recommendation = (
                        "Transfer from assigned intermediate."
                    )

            else:

                source_type = "Review"
                source_id = "No Valid Intermediate"
                source_conc_mM = None
                transfer_nL = direct_transfer_nL
                source_solvent_percent_for_assay = source_stock_solvent_percent
                assignment_status = "Review"

                recommendation = (
                    "Dose is not direct-transfer capable and no intermediate was assigned."
                )

            assay_backfill = calculate_solvent_backfill(
                transfer_nL=transfer_nL,
                final_assay_volume_uL=assay_well_volume_uL,
                target_solvent_percent=target_final_solvent_percent,
                source_solvent_percent=source_solvent_percent_for_assay
            )

            assignment_rows.append(
                {
                    "Compound ID": compound_id,
                    "Aliquot ID": aliquot_id,
                    "Dose Point": idx,
                    "Target Conc (mM)": target_mM,
                    "Target Conc": format_concentration(target_mM),
                    "Source Type": source_type,
                    "Source ID": source_id,
                    "Source Conc (mM)": source_conc_mM,
                    "Source Conc": format_concentration(source_conc_mM),
                    "Source Solvent %": source_solvent_percent_for_assay,
                    "Assay Transfer Numeric (nL)": transfer_nL,
                    "Assay Transfer": format_volume_nL(transfer_nL),
                    "Solvent From Transfer": format_volume_nL(
                        assay_backfill["Solvent From Transfer (nL)"]
                    ),
                    "Solvent Backfill": format_volume_nL(
                        assay_backfill["Solvent Backfill (nL)"]
                    ),
                    "Aqueous / Media Backfill": format_volume_nL(
                        assay_backfill["Aqueous Backfill (nL)"]
                    ),
                    "Actual Final Solvent %": assay_backfill[
                        "Actual Final Solvent %"
                    ],
                    "Assignment Status": assignment_status,
                    "Recommendation": recommendation
                }
            )

    intermediate_df = pd.DataFrame(
        intermediate_rows
    )

    assignment_df = pd.DataFrame(
        assignment_rows
    )

    return intermediate_df, assignment_df


# ==================================================
# PLATE MAP FUNCTIONS
# ==================================================

def get_plate_dimensions(plate_format):
    if plate_format == 96:
        rows = list(string.ascii_uppercase[:8])
        cols = list(range(1, 13))

    elif plate_format == 384:
        rows = list(string.ascii_uppercase[:16])
        cols = list(range(1, 25))

    elif plate_format == 1536:
        rows = [f"R{i + 1}" for i in range(32)]
        cols = list(range(1, 49))

    else:
        rows = list(string.ascii_uppercase[:16])
        cols = list(range(1, 25))

    return rows, cols


def normalize_well(well):
    well = str(well).strip().upper()

    match = re.match(
        r"^([A-Z]+)(\d+)$",
        well
    )

    if not match:
        return well

    row = match.group(1)
    col = match.group(2).zfill(2)

    return f"{row}{col}"


def parse_wells(text):
    if not text:
        return []

    return [
        normalize_well(x)
        for x in str(text).replace(";", ",").split(",")
        if x.strip()
    ]


def validate_wells(wells, plate_format):
    rows, cols = get_plate_dimensions(plate_format)
    col_labels = [str(c).zfill(2) for c in cols]

    valid = []
    invalid = []

    for well in wells:
        well = normalize_well(well)

        if len(well) < 3:
            invalid.append(well)
            continue

        row_part = well[:-2]
        col_part = well[-2:]

        if row_part in rows and col_part in col_labels:
            valid.append(well)
        else:
            invalid.append(well)

    return valid, invalid


def control_z_value(control_type):
    if control_type == "Vehicle Control":
        return 0.85

    if control_type == "Positive Control":
        return 0.92

    if control_type == "Negative Control":
        return 0.97

    if control_type == "Blank":
        return 1.0

    return 0.0


def build_control_rows(control_type, wells):
    control_rows = []

    for well in wells:
        well = normalize_well(well)
        row = well[:-2]
        col = well[-2:]

        control_rows.append(
            {
                "Well": well,
                "Row": row,
                "Column": col,
                "Well Type": control_type,
                "Compound ID": "",
                "Aliquot ID": "",
                "Dose Point": "",
                "Replicate": "",
                "Target Conc": "",
                "Source Type": control_type,
                "Source ID": control_type,
                "Source Conc": "",
                "Assay Transfer": "",
                "Solvent Backfill": "",
                "Aqueous / Media Backfill": "",
                "Actual Final Solvent %": "",
                "Assignment Status": "",
                "Map Label": control_type.replace(" Control", ""),
                "Hover Label": (
                    f"Well: {well}<br>"
                    f"Type: {control_type}"
                ),
                "Z Value": control_z_value(control_type)
            }
        )

    return control_rows


def build_automatic_control_wells(plate_format):
    rows, cols = get_plate_dimensions(plate_format)
    col_labels = [str(c).zfill(2) for c in cols]

    if plate_format == 96:
        control_row = rows[-1]

        vehicle = [
            f"{control_row}{col_labels[0]}",
            f"{control_row}{col_labels[1]}"
        ]

        positive = [
            f"{control_row}{col_labels[2]}",
            f"{control_row}{col_labels[3]}"
        ]

        negative = [
            f"{control_row}{col_labels[4]}",
            f"{control_row}{col_labels[5]}"
        ]

        blank = [
            f"{control_row}{col_labels[6]}",
            f"{control_row}{col_labels[7]}"
        ]

    elif plate_format == 384:
        vehicle = [f"{rows[-1]}23", f"{rows[-1]}24"]
        positive = [f"{rows[-2]}23", f"{rows[-2]}24"]
        negative = [f"{rows[-3]}23", f"{rows[-3]}24"]
        blank = [f"{rows[-4]}23", f"{rows[-4]}24"]

    else:
        last_col = col_labels[-1]
        second_last_col = col_labels[-2]

        vehicle = [
            f"{rows[-1]}{second_last_col}",
            f"{rows[-1]}{last_col}"
        ]

        positive = [
            f"{rows[-2]}{second_last_col}",
            f"{rows[-2]}{last_col}"
        ]

        negative = [
            f"{rows[-3]}{second_last_col}",
            f"{rows[-3]}{last_col}"
        ]

        blank = [
            f"{rows[-4]}{second_last_col}",
            f"{rows[-4]}{last_col}"
        ]

    return vehicle, positive, negative, blank


def build_dose_response_plate_map(
    dose_assignment_df,
    plate_format,
    replicates,
    start_column,
    reserved_wells
):
    rows, cols = get_plate_dimensions(plate_format)
    col_labels = [str(c).zfill(2) for c in cols]

    start_col_index = int(start_column) - 1

    records = []

    compounds = (
        dose_assignment_df["Compound ID"]
        .dropna()
        .unique()
        .tolist()
    )

    if len(dose_assignment_df) > 0:
        max_dose = int(dose_assignment_df["Dose Point"].max())
    else:
        max_dose = 1

    current_row_idx = 0

    for compound_id in compounds:

        compound_df = (
            dose_assignment_df[
                dose_assignment_df["Compound ID"] == compound_id
            ]
            .sort_values("Dose Point")
        )

        for rep in range(1, int(replicates) + 1):

            if current_row_idx >= len(rows):
                break

            row_label = rows[current_row_idx]

            for _, dose_row in compound_df.iterrows():

                dose_point = int(dose_row["Dose Point"])
                col_index = start_col_index + dose_point - 1

                if col_index >= len(col_labels):
                    continue

                col_label = col_labels[col_index]
                well = f"{row_label}{col_label}"

                if well in reserved_wells:
                    continue

                dose_intensity = (
                    (max_dose - dose_point + 1)
                    / max_dose
                )

                records.append(
                    {
                        "Well": well,
                        "Row": row_label,
                        "Column": col_label,
                        "Well Type": "Compound",
                        "Compound ID": compound_id,
                        "Aliquot ID": dose_row.get("Aliquot ID", ""),
                        "Dose Point": dose_point,
                        "Replicate": rep,
                        "Target Conc": dose_row.get("Target Conc", ""),
                        "Source Type": dose_row.get("Source Type", ""),
                        "Source ID": dose_row.get("Source ID", ""),
                        "Source Conc": dose_row.get("Source Conc", ""),
                        "Assay Transfer": dose_row.get("Assay Transfer", ""),
                        "Solvent Backfill": dose_row.get("Solvent Backfill", ""),
                        "Aqueous / Media Backfill": dose_row.get(
                            "Aqueous / Media Backfill",
                            ""
                        ),
                        "Actual Final Solvent %": dose_row.get(
                            "Actual Final Solvent %",
                            ""
                        ),
                        "Assignment Status": dose_row.get(
                            "Assignment Status",
                            ""
                        ),
                        "Map Label": f"D{dose_point}<br>R{rep}",
                        "Hover Label": (
                            f"Well: {well}<br>"
                            f"Type: Compound<br>"
                            f"Compound: {compound_id}<br>"
                            f"Aliquot: {dose_row.get('Aliquot ID', '')}<br>"
                            f"Dose Point: D{dose_point}<br>"
                            f"Replicate: {rep}<br>"
                            f"Target Conc: {dose_row.get('Target Conc', '')}<br>"
                            f"Source: {dose_row.get('Source Type', '')} - "
                            f"{dose_row.get('Source ID', '')}<br>"
                            f"Transfer: {dose_row.get('Assay Transfer', '')}<br>"
                            f"Solvent Backfill: {dose_row.get('Solvent Backfill', '')}<br>"
                            f"Media Backfill: {dose_row.get('Aqueous / Media Backfill', '')}<br>"
                            f"Assignment: {dose_row.get('Assignment Status', '')}"
                        ),
                        "Z Value": 0.10 + dose_intensity * 0.65
                    }
                )

            current_row_idx += 1

    return pd.DataFrame(records)


def add_control_wells_to_plate_map(
    plate_map_df,
    vehicle_wells,
    positive_wells,
    negative_wells,
    blank_wells
):
    control_rows = []

    control_rows.extend(
        build_control_rows(
            "Vehicle Control",
            vehicle_wells
        )
    )

    control_rows.extend(
        build_control_rows(
            "Positive Control",
            positive_wells
        )
    )

    control_rows.extend(
        build_control_rows(
            "Negative Control",
            negative_wells
        )
    )

    control_rows.extend(
        build_control_rows(
            "Blank",
            blank_wells
        )
    )

    if len(control_rows) == 0:
        return plate_map_df

    control_df = pd.DataFrame(control_rows)

    return pd.concat(
        [plate_map_df, control_df],
        ignore_index=True
    )


def empty_plate_grids(plate_format):
    rows, cols = get_plate_dimensions(plate_format)
    col_labels = [str(c).zfill(2) for c in cols]

    label_grid = pd.DataFrame(
        "",
        index=rows,
        columns=col_labels
    )

    hover_grid = pd.DataFrame(
        "",
        index=rows,
        columns=col_labels
    )

    z_grid = pd.DataFrame(
        0.0,
        index=rows,
        columns=col_labels
    )

    return label_grid, hover_grid, z_grid


def plot_dose_response_plate_map(
    plate_map_df,
    plate_format,
    title
):
    label_grid, hover_grid, z_grid = empty_plate_grids(
        plate_format
    )

    for _, row in plate_map_df.iterrows():
        r = row["Row"]
        c = row["Column"]

        if r in label_grid.index and c in label_grid.columns:
            label_grid.loc[r, c] = row["Map Label"]
            hover_grid.loc[r, c] = row["Hover Label"]
            z_grid.loc[r, c] = row["Z Value"]

    colorscale = [
        [0.00, "#F9FAFB"],
        [0.10, "#DBEAFE"],
        [0.35, "#60A5FA"],
        [0.70, "#1D4ED8"],
        [0.78, "#9CA3AF"],
        [0.86, "#22C55E"],
        [0.94, "#EF4444"],
        [1.00, "#FACC15"]
    ]

    fig = go.Figure(
        data=go.Heatmap(
            z=z_grid.values,
            x=z_grid.columns,
            y=z_grid.index,
            text=label_grid.values,
            hovertext=hover_grid.values,
            hoverinfo="text",
            colorscale=colorscale,
            showscale=False
        )
    )

    fig.update_traces(
        texttemplate="%{text}",
        textfont={
            "size": 8,
            "color": "black"
        }
    )

    fig.update_layout(
        title=title,
        xaxis_title="Column",
        yaxis_title="Row",
        height=760,
        plot_bgcolor="white"
    )

    fig.update_yaxes(
        autorange="reversed"
    )

    return fig


def convert_df_to_excel_bytes(df):
    output = BytesIO()

    with pd.ExcelWriter(
        output,
        engine="openpyxl"
    ) as writer:
        df.to_excel(
            writer,
            index=False,
            sheet_name="Dose Response Plate Map"
        )

    return output.getvalue()


# ==================================================
# TABS
# ==================================================

tab_inventory, tab_assay, tab_intermediate, tab_plate = st.tabs(
    [
        "Inventory Search",
        "Assay Designer",
        "Intermediate Strategy",
        "Dose Response Plate Map"
    ]
)


# ==================================================
# TAB 1: INVENTORY SEARCH WITH GLOBAL SEARCH + FILTERS
# ==================================================

with tab_inventory:

    st.subheader("Inventory Search")

    st.write(
        """
        Search inventory by compound ID, solvent, container type, status,
        location, or any other visible inventory field. Then select the exact
        aliquot/container record to send into assay design.
        """
    )

    search_term = st.text_input(
        "Global Search",
        placeholder="Examples: CMPD-1001, DMSO, Tube, Available, Freezer A, 10",
        key="global_inventory_search"
    )

    with st.expander("Optional Filters", expanded=True):

        f1, f2, f3, f4 = st.columns(4)

        with f1:
            solvent_filter = st.multiselect(
                "Solvent",
                clean_unique_values(inventory["solvent"]),
                key="solvent_filter"
            )

        with f2:
            container_filter = st.multiselect(
                "Container Type",
                clean_unique_values(inventory["container_type"]),
                key="container_filter"
            )

        with f3:
            status_filter = st.multiselect(
                "Status",
                clean_unique_values(inventory["inventory_status"]),
                key="status_filter"
            )

        with f4:
            location_filter = st.multiselect(
                "Location",
                clean_unique_values(inventory["location"]),
                key="location_filter"
            )

    run_search = st.button(
        "Search Inventory",
        key="inventory_search_button"
    )

    if run_search:

        search_results = run_inventory_global_search(
            inventory_df=inventory,
            search_term=search_term,
            solvent_filter=solvent_filter,
            container_filter=container_filter,
            status_filter=status_filter,
            location_filter=location_filter
        )

        if len(search_results) == 0:

            st.warning(
                "No matching inventory records found."
            )

        else:

            search_results["usable_uL"] = (
                search_results["available_uL"] - 10
            )

            search_results["Inventory_Health"] = (
                search_results.apply(
                    inventory_health,
                    axis=1
                )
            )

            search_results["SOP_Guidance"] = (
                search_results["solvent"]
                .apply(
                    solvent_guidance
                )
            )

            search_results["Use"] = False

            st.session_state[
                "last_inventory_search_results"
            ] = search_results

            st.success(
                f"{len(search_results)} matching inventory record(s) found."
            )

    if "last_inventory_search_results" in st.session_state:

        search_results = st.session_state[
            "last_inventory_search_results"
        ]

        st.subheader(
            "Select Aliquots / Inventory Records to Use"
        )

        selection_cols = [
            "Use",
            "Aliquot ID",
            "compound_id",
            "location",
            "container_type",
            "stock_mM",
            "available_uL",
            "usable_uL",
            "solvent",
            "freeze_thaw",
            "inventory_status",
            "expiry",
            "Inventory_Health"
        ]

        selection_df = st.data_editor(
            search_results[selection_cols],
            hide_index=True,
            use_container_width=True,
            column_config={
                "Use": st.column_config.CheckboxColumn(
                    "Use This Aliquot",
                    help="Select the aliquot/container record to send into assay design."
                ),
                "Aliquot ID": st.column_config.TextColumn(
                    "Aliquot ID",
                    disabled=True
                )
            },
            disabled=[
                "Aliquot ID",
                "compound_id",
                "location",
                "container_type",
                "stock_mM",
                "available_uL",
                "usable_uL",
                "solvent",
                "freeze_thaw",
                "inventory_status",
                "expiry",
                "Inventory_Health"
            ],
            key="aliquot_selection_editor"
        )

        proceed_selected = st.button(
            "Use Selected Aliquots",
            key="use_selected_aliquots_button"
        )

        if proceed_selected:

            selected_ids = (
                selection_df[
                    selection_df["Use"] == True
                ]["Aliquot ID"]
                .tolist()
            )

            selected_compounds = (
                search_results[
                    search_results["Aliquot ID"].isin(selected_ids)
                ]
                .copy()
            )

            if len(selected_compounds) == 0:

                st.warning(
                    "Select at least one aliquot or inventory record."
                )

            else:

                selected_compounds["usable_uL"] = (
                    selected_compounds["available_uL"] - 10
                )

                selected_compounds["Inventory_Health"] = (
                    selected_compounds.apply(
                        inventory_health,
                        axis=1
                    )
                )

                selected_compounds["SOP_Guidance"] = (
                    selected_compounds["solvent"]
                    .apply(
                        solvent_guidance
                    )
                )

                st.session_state[
                    "selected_compounds"
                ] = selected_compounds

                st.success(
                    f"{len(selected_compounds)} selected aliquot(s) loaded for assay design."
                )

        if "selected_compounds" in st.session_state:

            st.subheader(
                "Selected Aliquots for Assay Design"
            )

            selected_view_cols = [
                "Aliquot ID",
                "compound_id",
                "location",
                "container_type",
                "stock_mM",
                "available_uL",
                "usable_uL",
                "solvent",
                "freeze_thaw",
                "inventory_status",
                "expiry",
                "Inventory_Health"
            ]

            st.dataframe(
                st.session_state["selected_compounds"][selected_view_cols],
                use_container_width=True
            )

            st.subheader(
                "SOP Guidance for Selected Aliquots"
            )

            for _, row in st.session_state["selected_compounds"].iterrows():

                with st.expander(
                    f"{row['compound_id']} | {row['location']} | {row['container_type']}"
                ):

                    st.markdown(
                        f"**Aliquot ID:** {row['Aliquot ID']}"
                    )

                    st.markdown(
                        f"**Inventory Health:** {row['Inventory_Health']}"
                    )

                    st.markdown(
                        row["SOP_Guidance"]
                    )

    with st.expander(
        "View Full Inventory"
    ):

        st.dataframe(
            inventory,
            use_container_width=True
        )


# ==================================================
# TAB 2: ASSAY DESIGNER
# ==================================================

with tab_assay:

    st.subheader("Assay Designer / Direct Stock Calculator")

    st.write(
        """
        Build a dose-response curve, calculate direct-from-stock transfers,
        and calculate matching solvent backfill so every dose point reaches
        the same target final solvent percentage.
        """
    )

    if "selected_compounds" not in st.session_state:

        st.info(
            "Search inventory and select one or more aliquots first."
        )

    else:

        selected_compounds = st.session_state[
            "selected_compounds"
        ]

        st.subheader(
            "Selected Aliquots"
        )

        st.dataframe(
            selected_compounds[
                [
                    "Aliquot ID",
                    "compound_id",
                    "stock_mM",
                    "available_uL",
                    "usable_uL",
                    "solvent",
                    "inventory_status",
                    "freeze_thaw",
                    "expiry",
                    "location",
                    "container_type"
                ]
            ],
            use_container_width=True
        )

        st.divider()

        st.subheader(
            "Assay Parameters"
        )

        a1, a2, a3, a4 = st.columns(4)

        with a1:
            assay_type = st.selectbox(
                "Assay Type",
                [
                    "Dose Response",
                    "Single Point"
                ],
                key="assay_type"
            )

        with a2:
            top_conc_value = st.number_input(
                "Top Concentration",
                value=10.0,
                min_value=0.0,
                key="top_concentration"
            )

        with a3:
            top_conc_unit = st.selectbox(
                "Top Concentration Unit",
                [
                    "µM",
                    "mM",
                    "nM",
                    "pM"
                ],
                key="top_concentration_unit"
            )

        with a4:
            dilution_factor = st.number_input(
                "Dilution Factor",
                value=3.0,
                min_value=1.01,
                key="dilution_factor"
            )

        b1, b2, b3, b4 = st.columns(4)

        with b1:
            if assay_type == "Single Point":
                dose_points = 1

                st.number_input(
                    "Dose Points",
                    value=1,
                    disabled=True,
                    key="dose_points_disabled"
                )

            else:
                dose_points = st.number_input(
                    "Dose Points",
                    value=10,
                    min_value=1,
                    step=1,
                    key="dose_points"
                )

        with b2:
            assay_well_volume_uL = st.number_input(
                "Final Assay Well Volume (µL)",
                value=10.0,
                min_value=0.1,
                key="assay_well_volume"
            )

        with b3:
            plates = st.number_input(
                "Number of Assay Plates",
                value=1,
                min_value=1,
                step=1,
                key="assay_plates"
            )

        with b4:
            replicates = st.number_input(
                "Replicates",
                value=1,
                min_value=1,
                step=1,
                key="assay_replicates"
            )

        c1, c2, c3, c4 = st.columns(4)

        with c1:
            min_transfer_nL = st.number_input(
                "Minimum Echo Transfer (nL)",
                value=2.5,
                min_value=0.0,
                key="min_echo_transfer"
            )

        with c2:
            max_transfer_nL = st.number_input(
                "Maximum Echo Transfer (nL)",
                value=500.0,
                min_value=0.1,
                key="max_echo_transfer"
            )

        with c3:
            target_solvent_percent = st.number_input(
                "Target Final Solvent (%)",
                value=0.5,
                min_value=0.0,
                key="target_final_solvent"
            )

        with c4:
            source_stock_solvent_percent = st.number_input(
                "Source Stock Solvent (%)",
                value=100.0,
                min_value=0.0,
                max_value=100.0,
                key="source_stock_solvent_percent"
            )

        run_assay = st.button(
            "Generate Assay Design",
            key="generate_assay_design"
        )

        if run_assay:

            top_conc_mM = convert_concentration_to_mM(
                top_conc_value,
                top_conc_unit
            )

            concentrations_mM = build_concentration_series(
                top_conc_mM=top_conc_mM,
                dose_points=dose_points,
                dilution_factor=dilution_factor
            )

            concentration_df = pd.DataFrame(
                {
                    "Dose Point": range(
                        1,
                        int(dose_points) + 1
                    ),
                    "Concentration (mM)": concentrations_mM
                }
            )

            concentration_df["Display"] = (
                concentration_df["Concentration (mM)"]
                .apply(format_concentration)
            )

            st.session_state[
                "concentration_df"
            ] = concentration_df

            assay_params = {
                "assay_well_volume_uL": assay_well_volume_uL,
                "min_transfer_nL": min_transfer_nL,
                "max_transfer_nL": max_transfer_nL,
                "plates": plates,
                "replicates": replicates,
                "target_solvent_percent": target_solvent_percent,
                "source_stock_solvent_percent": source_stock_solvent_percent
            }

            st.session_state[
                "assay_params"
            ] = assay_params

            st.subheader(
                "Dose Response Curve"
            )

            fig = px.line(
                concentration_df,
                x="Dose Point",
                y="Concentration (mM)",
                markers=True,
                title="Assay Concentration Series"
            )

            fig.update_yaxes(
                type="log",
                title_text="Concentration (mM, log scale)"
            )

            fig.update_xaxes(
                title_text="Dose Point"
            )

            fig.update_layout(
                height=500
            )

            st.plotly_chart(
                fig,
                use_container_width=True
            )

            st.subheader(
                "Concentration Ladder"
            )

            st.dataframe(
                concentration_df[
                    [
                        "Dose Point",
                        "Display",
                        "Concentration (mM)"
                    ]
                ],
                use_container_width=True
            )

            direct_df = build_direct_stock_assay_design(
                selected_compounds=selected_compounds,
                concentrations_mM=concentrations_mM,
                assay_well_volume_uL=assay_well_volume_uL,
                min_transfer_nL=min_transfer_nL,
                max_transfer_nL=max_transfer_nL,
                plates=plates,
                replicates=replicates,
                target_solvent_percent=target_solvent_percent,
                source_solvent_percent=source_stock_solvent_percent
            )

            st.session_state[
                "direct_assay_df"
            ] = direct_df

            st.subheader(
                "Direct Stock Transfer and Solvent Backfill Review"
            )

            review_count = len(
                direct_df[
                    direct_df["Assay Status"] == "Review"
                ]
            )

            ok_count = len(
                direct_df[
                    direct_df["Assay Status"] == "Direct OK"
                ]
            )

            s1, s2, s3 = st.columns(3)

            s1.metric(
                "Dose Rows",
                len(direct_df)
            )

            s2.metric(
                "Direct OK",
                ok_count
            )

            s3.metric(
                "Review",
                review_count
            )

            display_cols = [
                "Compound ID",
                "Aliquot ID",
                "Dose Point",
                "Target Conc",
                "Stock (mM)",
                "Solvent",
                "Source Solvent %",
                "Compound Transfer",
                "Target Total Solvent",
                "Solvent From Compound",
                "Solvent Backfill",
                "Aqueous / Media Backfill",
                "Actual Final Solvent %",
                "Required Source Volume",
                "Usable Inventory (uL)",
                "Assay Status",
                "Review Reasons",
                "Recommendation"
            ]

            st.dataframe(
                direct_df[display_cols],
                use_container_width=True
            )


# ==================================================
# TAB 3: INTERMEDIATE STRATEGY
# ==================================================

with tab_intermediate:

    st.subheader("Intermediate Strategy Designer")

    if (
        "selected_compounds" not in st.session_state
        or "concentration_df" not in st.session_state
        or "assay_params" not in st.session_state
    ):

        st.info(
            "Complete Inventory Search and Assay Designer first."
        )

    else:

        selected_compounds = st.session_state[
            "selected_compounds"
        ]

        concentration_df = st.session_state[
            "concentration_df"
        ]

        assay_params = st.session_state[
            "assay_params"
        ]

        concentrations_mM = concentration_df[
            "Concentration (mM)"
        ].tolist()

        left, right = st.columns(2)

        with left:
            st.write("Selected Aliquots")

            st.dataframe(
                selected_compounds[
                    [
                        "Aliquot ID",
                        "compound_id",
                        "stock_mM",
                        "available_uL",
                        "solvent",
                        "location"
                    ]
                ],
                use_container_width=True
            )

        with right:
            st.write("Concentration Ladder")

            st.dataframe(
                concentration_df[
                    [
                        "Dose Point",
                        "Display",
                        "Concentration (mM)"
                    ]
                ],
                use_container_width=True
            )

        st.divider()

        st.subheader(
            "Intermediate Strategy Parameters"
        )

        m1, m2, m3 = st.columns(3)

        with m1:
            intermediate_start_mode = st.selectbox(
                "Intermediate Start Mode",
                [
                    "Auto from first non-direct dose",
                    "Manual start dose point"
                ],
                key="intermediate_start_mode"
            )

            max_dose_point = max(
                1,
                len(concentrations_mM)
            )

            default_start_dose = min(
                4,
                max_dose_point
            )

            manual_start_dose = st.number_input(
                "Manual Start Dose Point",
                value=default_start_dose,
                min_value=1,
                max_value=max_dose_point,
                step=1,
                key="manual_start_dose"
            )

        with m2:
            max_intermediates_per_compound = st.number_input(
                "Max Intermediates Per Compound",
                value=3,
                min_value=1,
                step=1,
                key="max_intermediates"
            )

            preferred_intermediate_transfer_nL = st.number_input(
                "Preferred Intermediate Transfer (nL)",
                value=10.0,
                min_value=0.1,
                key="preferred_intermediate_transfer"
            )

        with m3:
            intermediate_working_volume_uL = st.number_input(
                "Intermediate Working Volume (µL)",
                value=5.0,
                min_value=0.1,
                key="intermediate_working_volume"
            )

            intermediate_dead_volume_uL = st.number_input(
                "Intermediate Dead Volume (µL)",
                value=2.5,
                min_value=0.0,
                key="intermediate_dead_volume"
            )

        n1, n2 = st.columns(2)

        with n1:
            intermediate_solvent_percent = st.number_input(
                "Target Intermediate Solvent (%)",
                value=10.0,
                min_value=0.0,
                max_value=100.0,
                key="intermediate_solvent_percent"
            )

        with n2:
            dilution_buffer_name = st.text_input(
                "Diluent / Buffer",
                value="Assay Buffer",
                key="dilution_buffer_name"
            )

        run_intermediate = st.button(
            "Generate Intermediate Strategy",
            key="generate_intermediate_strategy"
        )

        if run_intermediate:

            intermediate_df, assignment_df = build_intermediate_strategy(
                selected_compounds=selected_compounds,
                concentrations_mM=concentrations_mM,
                assay_well_volume_uL=assay_params["assay_well_volume_uL"],
                min_transfer_nL=assay_params["min_transfer_nL"],
                max_transfer_nL=assay_params["max_transfer_nL"],
                target_final_solvent_percent=assay_params["target_solvent_percent"],
                source_stock_solvent_percent=assay_params["source_stock_solvent_percent"],
                preferred_intermediate_transfer_nL=preferred_intermediate_transfer_nL,
                max_intermediates_per_compound=max_intermediates_per_compound,
                intermediate_working_volume_uL=intermediate_working_volume_uL,
                intermediate_dead_volume_uL=intermediate_dead_volume_uL,
                intermediate_solvent_percent=intermediate_solvent_percent,
                dilution_buffer_name=dilution_buffer_name,
                intermediate_start_mode=intermediate_start_mode,
                manual_start_dose=manual_start_dose
            )

            st.session_state[
                "intermediate_df"
            ] = intermediate_df

            st.session_state[
                "dose_assignment_df"
            ] = assignment_df

            st.success(
                "Intermediate strategy generated."
            )

            st.subheader(
                "Intermediate Recipes"
            )

            if len(intermediate_df) == 0:

                st.info(
                    "No intermediates were generated."
                )

            else:

                recipe_cols = [
                    "Compound ID",
                    "Intermediate ID",
                    "Intermediate Conc",
                    "Source Stock (mM)",
                    "Solvent",
                    "Intermediate Target Solvent %",
                    "Actual Intermediate Solvent %",
                    "Working Volume (uL)",
                    "Dead Volume (uL)",
                    "Total Intermediate Volume (uL)",
                    "Source Stock to Make Intermediate (uL)",
                    "Matching Solvent Backfill to Make Intermediate (uL)",
                    "Diluent / Buffer",
                    "Diluent Volume to Make Intermediate (uL)",
                    "Dose Points Covered",
                    "Recipe Flags",
                    "Recipe Recommendation"
                ]

                st.dataframe(
                    intermediate_df[
                        recipe_cols
                    ],
                    use_container_width=True
                )

            st.subheader(
                "Dose Source Assignment and Final Assay Backfill"
            )

            assignment_cols = [
                "Compound ID",
                "Aliquot ID",
                "Dose Point",
                "Target Conc",
                "Source Type",
                "Source ID",
                "Source Conc",
                "Source Solvent %",
                "Assay Transfer",
                "Solvent From Transfer",
                "Solvent Backfill",
                "Aqueous / Media Backfill",
                "Actual Final Solvent %",
                "Assignment Status",
                "Recommendation"
            ]

            st.dataframe(
                assignment_df[
                    assignment_cols
                ],
                use_container_width=True
            )

            st.subheader(
                "Intermediate Summary"
            )

            if len(intermediate_df) > 0:

                summary = (
                    intermediate_df
                    .groupby("Compound ID")
                    .agg(
                        Intermediate_Count=(
                            "Intermediate ID",
                            "nunique"
                        ),
                        Total_Source_Stock_uL=(
                            "Source Stock to Make Intermediate (uL)",
                            "sum"
                        ),
                        Total_Solvent_Backfill_uL=(
                            "Matching Solvent Backfill to Make Intermediate (uL)",
                            "sum"
                        ),
                        Total_Diluent_uL=(
                            "Diluent Volume to Make Intermediate (uL)",
                            "sum"
                        )
                    )
                    .reset_index()
                )

                st.dataframe(
                    summary,
                    use_container_width=True
                )


# ==================================================
# TAB 4: DOSE RESPONSE PLATE MAP
# ==================================================

with tab_plate:

    st.subheader("Dose Response Plate Map")

    if "dose_assignment_df" not in st.session_state:

        st.info(
            "Complete Intermediate Strategy first."
        )

    else:

        dose_assignment_df = st.session_state[
            "dose_assignment_df"
        ]

        p1, p2, p3 = st.columns(3)

        with p1:
            plate_format = st.selectbox(
                "Plate Format",
                [
                    96,
                    384,
                    1536
                ],
                index=1,
                key="plate_format"
            )

        with p2:
            map_replicates = st.number_input(
                "Replicates to Map",
                value=3,
                min_value=1,
                step=1,
                key="plate_map_replicates"
            )

        with p3:
            start_column = st.number_input(
                "Start Column",
                value=1,
                min_value=1,
                step=1,
                key="plate_start_column"
            )

        st.divider()

        st.subheader("Control Well Placement")

        control_mode = st.radio(
            "Control Placement Mode",
            [
                "Automatic Controls",
                "Manual Well Selection",
                "No Controls"
            ],
            key="control_placement_mode"
        )

        vehicle_wells = ""
        positive_wells = ""
        negative_wells = ""
        blank_wells = ""

        if control_mode == "Manual Well Selection":

            c1, c2 = st.columns(2)

            with c1:
                vehicle_wells = st.text_input(
                    "Vehicle Control Wells",
                    value="P23,P24",
                    key="vehicle_control_wells"
                )

                positive_wells = st.text_input(
                    "Positive Control Wells",
                    value="O23,O24",
                    key="positive_control_wells"
                )

            with c2:
                negative_wells = st.text_input(
                    "Negative Control Wells",
                    value="N23,N24",
                    key="negative_control_wells"
                )

                blank_wells = st.text_input(
                    "Blank Wells",
                    value="M23,M24",
                    key="blank_control_wells"
                )

        generate_map = st.button(
            "Generate Dose Response Plate Map",
            key="generate_plate_map"
        )

        if generate_map:

            if control_mode == "Manual Well Selection":

                vehicle_raw = parse_wells(
                    vehicle_wells
                )

                positive_raw = parse_wells(
                    positive_wells
                )

                negative_raw = parse_wells(
                    negative_wells
                )

                blank_raw = parse_wells(
                    blank_wells
                )

                vehicle_valid, vehicle_invalid = validate_wells(
                    vehicle_raw,
                    plate_format
                )

                positive_valid, positive_invalid = validate_wells(
                    positive_raw,
                    plate_format
                )

                negative_valid, negative_invalid = validate_wells(
                    negative_raw,
                    plate_format
                )

                blank_valid, blank_invalid = validate_wells(
                    blank_raw,
                    plate_format
                )

                invalid_all = (
                    vehicle_invalid
                    + positive_invalid
                    + negative_invalid
                    + blank_invalid
                )

                if invalid_all:
                    st.warning(
                        "These wells were not valid for the selected plate format: "
                        + ", ".join(invalid_all)
                    )

                vehicle_set = set(vehicle_valid)
                positive_set = set(positive_valid)
                negative_set = set(negative_valid)
                blank_set = set(blank_valid)

            elif control_mode == "Automatic Controls":

                vehicle_auto, positive_auto, negative_auto, blank_auto = (
                    build_automatic_control_wells(
                        plate_format
                    )
                )

                vehicle_set = set(vehicle_auto)
                positive_set = set(positive_auto)
                negative_set = set(negative_auto)
                blank_set = set(blank_auto)

            else:

                vehicle_set = set()
                positive_set = set()
                negative_set = set()
                blank_set = set()

            reserved_wells = (
                vehicle_set
                | positive_set
                | negative_set
                | blank_set
            )

            plate_map_df = build_dose_response_plate_map(
                dose_assignment_df=dose_assignment_df,
                plate_format=plate_format,
                replicates=map_replicates,
                start_column=start_column,
                reserved_wells=reserved_wells
            )

            plate_map_df = add_control_wells_to_plate_map(
                plate_map_df=plate_map_df,
                vehicle_wells=sorted(vehicle_set),
                positive_wells=sorted(positive_set),
                negative_wells=sorted(negative_set),
                blank_wells=sorted(blank_set)
            )

            st.session_state[
                "plate_map_df"
            ] = plate_map_df

            total_wells = (
                len(get_plate_dimensions(plate_format)[0])
                * len(get_plate_dimensions(plate_format)[1])
            )

            used_wells = len(plate_map_df)
            remaining_wells = total_wells - used_wells

            st.subheader("Plate Summary")

            s1, s2, s3, s4 = st.columns(4)

            s1.metric(
                "Plate Format",
                f"{plate_format}-well"
            )

            s2.metric(
                "Used Wells",
                used_wells
            )

            s3.metric(
                "Remaining Wells",
                remaining_wells
            )

            s4.metric(
                "Unique Compounds",
                plate_map_df[
                    plate_map_df["Well Type"] == "Compound"
                ]["Compound ID"].nunique()
            )

            if remaining_wells < 0:
                st.error(
                    "Plate map exceeds available wells. Reduce compounds, dose points, replicates, or controls."
                )

            else:

                st.subheader("Visual Dose Response Plate Map")

                fig = plot_dose_response_plate_map(
                    plate_map_df=plate_map_df,
                    plate_format=plate_format,
                    title=f"{plate_format}-Well Dose Response Plate Map"
                )

                st.plotly_chart(
                    fig,
                    use_container_width=True
                )

                st.subheader("Plate Map Table")

                table_cols = [
                    "Well",
                    "Well Type",
                    "Compound ID",
                    "Aliquot ID",
                    "Dose Point",
                    "Replicate",
                    "Target Conc",
                    "Source Type",
                    "Source ID",
                    "Assay Transfer",
                    "Solvent Backfill",
                    "Aqueous / Media Backfill",
                    "Actual Final Solvent %",
                    "Assignment Status"
                ]

                st.dataframe(
                    plate_map_df[table_cols],
                    use_container_width=True
                )

                st.subheader("Export Plate Map")

                csv_bytes = (
                    plate_map_df
                    .to_csv(index=False)
                    .encode("utf-8")
                )

                st.download_button(
                    label="Download Plate Map CSV",
                    data=csv_bytes,
                    file_name="dose_response_plate_map.csv",
                    mime="text/csv"
                )

                xlsx_bytes = convert_df_to_excel_bytes(
                    plate_map_df
                )

                st.download_button(
                    label="Download Plate Map XLSX",
                    data=xlsx_bytes,
                    file_name="dose_response_plate_map.xlsx",
                    mime=(
                        "application/vnd.openxmlformats-officedocument."
                        "spreadsheetml.sheet"
                    )
                )