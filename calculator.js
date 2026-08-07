/* ============================================================
   PEPTIDE FRAGMENT CALCULATOR for Mass Spectrometry
   ============================================================ */

"use strict";

// ─── Constants ───────────────────────────────────────────────
const PROTON = 1.007276;
const WATER_MONO = 18.010565;
const WATER_AVG  = 18.01528;
const CO   = 27.994915;
const NH3  = 17.026549;

const MONO = {
  G:57.02146,  A:71.03711,  V:99.06841,  L:113.08406, I:113.08406,
  P:97.05276,  F:147.06841, W:186.07931, M:131.04049, S:87.03203,
  T:101.04768, C:103.00919, Y:163.06333, H:137.05891, D:115.02694,
  E:129.04259, N:114.04293, Q:128.05858, K:128.09496, R:156.10111
};

const AVG = {
  G:57.0519,  A:71.0788,  V:99.1326,  L:113.1594, I:113.1594,
  P:97.1167,  F:147.1766, W:186.2132, M:131.1926, S:87.0782,
  T:101.1051, C:103.1388, Y:163.1760, H:137.1411, D:115.0886,
  E:129.1155, N:114.1038, Q:128.1307, K:128.1741, R:156.1875
};

// Static modification definitions
// target: single AA letter | 'nterm' | 'cterm'
const BUILTIN_MODS = [
  { name: "Acetylation of K",          target: "K",     mass: 42.0106  },
  { name: "Acetylation of nterm",       target: "nterm", mass: 42.0106  },
  { name: "Carbamidomethyl C",          target: "C",     mass: 57.0215  },
  { name: "Cysteinylated C",            target: "C",     mass: 119.0041 },
  { name: "Dimethylation of K",         target: "K",     mass: 28.0313  },
  { name: "Dimethylation of R",         target: "R",     mass: 28.0313  },
  { name: "GlcNac on S",               target: "S",     mass: 203.0794 },
  { name: "GlcNac on T",               target: "T",     mass: 203.0794 },
  { name: "GlcNac+Galnaz+PCBiotin on S",target: "S",    mass: 901.3416 },
  { name: "Heavy R",                   target: "R",     mass: 10.0083  },
  { name: "Methylation of cterm",       target: "cterm", mass: 14.0157  },
  { name: "Methylation of D",          target: "D",     mass: 14.0157  },
  { name: "Methylation of E",          target: "E",     mass: 14.0157  },
  { name: "Methylation of K",          target: "K",     mass: 14.0157  },
  { name: "Methylation of nterm",      target: "nterm", mass: 14.0157  },
  { name: "Methylation of R",          target: "R",     mass: 14.0157  },
  { name: "Oxidation of M",            target: "M",     mass: 15.9949  },
  { name: "Phosphorylation of S",      target: "S",     mass: 79.9663  },
  { name: "Phosphorylation of T",      target: "T",     mass: 79.9663  },
  { name: "Phosphorylation of Y",      target: "Y",     mass: 79.9663  },
  { name: "Propionyl K",               target: "K",     mass: 56.0262  },
  { name: "Trimethylation of K",       target: "K",     mass: 42.0470  },
  { name: "Trimethylation of R",       target: "R",     mass: 42.0470  },
  { name: "C-terminus OH to NH2",       target: "cterm", mass: -0.9840  }
];

// All mods (built-in + user-added)
let allMods = [...BUILTIN_MODS];

// ─── DOM References ──────────────────────────────────────────
const peptideInput      = document.getElementById("peptide-input");
const precursorMinEl    = document.getElementById("precursor-min");
const precursorMaxEl    = document.getElementById("precursor-max");
const fragmentMinEl     = document.getElementById("fragment-min");
const fragmentMaxEl     = document.getElementById("fragment-max");
const ionB              = document.getElementById("ion-b");
const ionY              = document.getElementById("ion-y");
const ionA              = document.getElementById("ion-a");
const ionC              = document.getElementById("ion-c");
const ionZ              = document.getElementById("ion-z");
const staticModsEl      = document.getElementById("static-mods-container");
const customNameEl      = document.getElementById("custom-mod-name");
const customTargetEl    = document.getElementById("custom-mod-target");
const customMassEl      = document.getElementById("custom-mod-mass");
const addCustomBtn      = document.getElementById("add-custom-mod-btn");
const calculateBtn      = document.getElementById("calculate-btn");
const clearBtn          = document.getElementById("clear-btn");
const resultsSection    = document.getElementById("results-section");
const precursorSection  = document.getElementById("precursor-section");
const fragmentSection   = document.getElementById("fragment-section");
const precursorSeqEl    = document.getElementById("precursor-sequence-display");
const fragmentSeqEl     = document.getElementById("fragment-sequence-display");
const precursorTableEl  = document.getElementById("precursor-table-container");
const fragmentTableEl   = document.getElementById("fragment-table-container");
const pageBreakEl       = document.getElementById("page-break");
const exportPdfBtn      = document.getElementById("export-pdf-btn");
const exportCsvBtn      = document.getElementById("export-csv-btn");

// ─── Initialise static mods list ────────────────────────────
function renderModsList() {
  staticModsEl.innerHTML = "";
  allMods.forEach((mod, idx) => {
    const div = document.createElement("div");
    div.className = "mod-item";
    div.innerHTML =
      `<input type="checkbox" id="mod-${idx}" data-idx="${idx}"/>` +
      `<label for="mod-${idx}">${mod.name} <span style="color:#94a3b8">(${mod.mass > 0 ? "+" : ""}${mod.mass})</span></label>`;
    staticModsEl.appendChild(div);
  });
}
renderModsList();

// ─── Add custom mod ──────────────────────────────────────────
addCustomBtn.addEventListener("click", () => {
  const name   = customNameEl.value.trim();
  const target = customTargetEl.value.trim().toUpperCase();
  const mass   = parseFloat(customMassEl.value);

  if (!name) { alert("Please enter a modification name."); return; }
  if (!target) { alert("Please enter a target (AA letter, nterm, or cterm)."); return; }
  if (isNaN(mass)) { alert("Please enter a valid mass shift."); return; }

  const normTarget = (target === "NTERM") ? "nterm" : (target === "CTERM") ? "cterm" : target[0];
  allMods.push({ name, target: normTarget, mass });
  renderModsList();
  customNameEl.value = "";
  customTargetEl.value = "";
  customMassEl.value = "";
});

// ─── Parse peptide sequence ───────────────────────────────────
// Returns array of { aa, monoMass, avgMass, label }
// Supports inline mods: K[+42.0106] or M[-16]
function parsePeptide(raw) {
  const seq = raw.trim().toUpperCase();
  const residues = [];
  let i = 0;
  while (i < seq.length) {
    const aa = seq[i];
    if (!MONO[aa]) {
      if (aa === "[") { i++; continue; } // skip stray bracket
      alert(`Unknown amino acid: "${aa}" at position ${i + 1}`);
      return null;
    }
    let inlineMod = 0;
    let displayMod = "";
    // Check for inline modification [+/-mass]
    if (i + 1 < seq.length && seq[i + 1] === "[") {
      const closeIdx = seq.indexOf("]", i + 2);
      if (closeIdx !== -1) {
        const modStr = seq.slice(i + 2, closeIdx);
        const parsed = parseFloat(modStr);
        if (!isNaN(parsed)) {
          inlineMod = parsed;
          displayMod = `[${modStr}]`;
          i = closeIdx + 1;
        } else { i++; }
      } else { i++; }
    } else { i++; }

    residues.push({
      aa,
      monoMass: MONO[aa] + inlineMod,
      avgMass:  AVG[aa]  + inlineMod,
      label: aa + displayMod
    });
  }
  return residues;
}

// ─── Apply static modifications ──────────────────────────────
function applyStaticMods(residues) {
  const checked = [...document.querySelectorAll("#static-mods-container input[type=checkbox]:checked")]
    .map(cb => allMods[parseInt(cb.dataset.idx)]);

  // nterm / cterm mods
  let ntermDelta = 0;
  let ctermDelta = 0;
  checked.forEach(mod => {
    if (mod.target === "nterm") ntermDelta += mod.mass;
    if (mod.target === "cterm") ctermDelta += mod.mass;
  });

  const modified = residues.map((r, idx) => {
    let extraMono = 0;
    let extraAvg  = 0;
    // nterm mod goes onto first residue
    if (idx === 0) { extraMono += ntermDelta; extraAvg += ntermDelta; }
    // cterm mod goes onto last residue
    if (idx === residues.length - 1) { extraMono += ctermDelta; extraAvg += ctermDelta; }
    // residue-specific mods
    checked.forEach(mod => {
      if (mod.target === r.aa) { extraMono += mod.mass; extraAvg += mod.mass; }
    });
    return {
      ...r,
      monoMass: r.monoMass + extraMono,
      avgMass:  r.avgMass  + extraAvg
    };
  });
  return modified;
}

// ─── Mass Calculation helpers ─────────────────────────────────
function precursorMono(residues, charge) {
  const pepMass = residues.reduce((s, r) => s + r.monoMass, 0) + WATER_MONO;
  return (pepMass + charge * PROTON) / charge;
}

function precursorAvg(residues, charge) {
  const pepMass = residues.reduce((s, r) => s + r.avgMass, 0) + WATER_AVG;
  return (pepMass + charge * PROTON) / charge;
}

// Cumulative mass arrays (N-terminal prefix sums, no water)
function nTermSumsMono(residues) {
  const s = [0];
  residues.forEach(r => s.push(s[s.length - 1] + r.monoMass));
  return s;
}
function nTermSumsAvg(residues) {
  const s = [0];
  residues.forEach(r => s.push(s[s.length - 1] + r.avgMass));
  return s;
}

// Ion m/z
// useAvg: only when charge >= 2 AND user selected average
function ionMz(neutralMass, charge, useAvgMode) {
  // useAvgMode already baked into neutralMass; just add protons
  return (neutralMass + charge * PROTON) / charge;
}

function fmt(val) {
  if (val === null || val === undefined || isNaN(val)) return "-";
  return val.toFixed(4);
}

// ─── Calculate ───────────────────────────────────────────────
calculateBtn.addEventListener("click", calculate);

function getMassType() {
  return document.querySelector('input[name="mass-type"]:checked').value; // 'mono' | 'average'
}

function calculate() {
  const rawSeq = peptideInput.value;
  if (!rawSeq.trim()) { alert("Please enter a peptide sequence."); return; }

  const residues = parsePeptide(rawSeq);
  if (!residues || residues.length === 0) return;

  const modResidues = applyStaticMods(residues);
  const N = modResidues.length;

  const pMin = parseInt(precursorMinEl.value) || 1;
  const pMax = parseInt(precursorMaxEl.value) || 3;
  const fMin = parseInt(fragmentMinEl.value)  || 1;
  const fMax = parseInt(fragmentMaxEl.value)  || 2;
  const massType = getMassType();

  const useB = ionB.checked;
  const useY = ionY.checked;
  const useA = ionA.checked;
  const useC = ionC.checked;
  const useZ = ionZ.checked;

  // ── Sequence display ──────────────────────────────────────
  const seqStr = modResidues.map(r => r.label).join("");
  precursorSeqEl.textContent = seqStr;
  fragmentSeqEl.textContent  = seqStr;

  // ── Precursor Table ───────────────────────────────────────
  precursorTableEl.innerHTML = buildPrecursorTable(modResidues, pMin, pMax);

  // ── Fragment Table ────────────────────────────────────────
  fragmentTableEl.innerHTML = buildFragmentTable(
    modResidues, N, fMin, fMax, massType, useB, useY, useA, useC, useZ
  );

  // ── Show results ──────────────────────────────────────────
  resultsSection.classList.remove("hidden");

  // Display mode
  const displayMode = document.querySelector('input[name="display-mode"]:checked').value;
  applyDisplayMode(displayMode);

  resultsSection.scrollIntoView({ behavior: "smooth" });
}

// ─── Display mode radios ──────────────────────────────────────
document.querySelectorAll('input[name="display-mode"]').forEach(r => {
  r.addEventListener("change", () => {
    if (!resultsSection.classList.contains("hidden")) {
      applyDisplayMode(r.value);
    }
  });
});

function applyDisplayMode(mode) {
  if (mode === "separate") {
    pageBreakEl.classList.remove("hidden");
    pageBreakEl.classList.add("page-break-active");
    fragmentSection.classList.add("page-break-print");
  } else {
    pageBreakEl.classList.add("hidden");
    pageBreakEl.classList.remove("page-break-active");
    fragmentSection.classList.remove("page-break-print");
  }
}

// ─── Build Precursor Table ────────────────────────────────────
function buildPrecursorTable(residues, pMin, pMax) {
  let html = `<table><thead><tr>
    <th class="center-col">Charge State</th>
    <th class="b-ion">Monoisotopic m/z (Da)</th>
    <th class="y-ion">Average m/z (Da)</th>
  </tr></thead><tbody>`;

  for (let z = pMin; z <= pMax; z++) {
    const mono = precursorMono(residues, z);
    const avg  = precursorAvg(residues, z);
    html += `<tr>
      <td class="idx-cell">+${z}</td>
      <td class="b-ion">${fmt(mono)}</td>
      <td class="y-ion">${fmt(avg)}</td>
    </tr>`;
  }
  html += `</tbody></table>`;
  return html;
}

// ─── Build Fragment Table ─────────────────────────────────────
function buildFragmentTable(residues, N, fMin, fMax, massType, useB, useY, useA, useC, useZ) {
  // Precompute prefix/suffix sums
  const monoPrefix = nTermSumsMono(residues); // monoPrefix[i] = sum of first i residues
  const avgPrefix  = nTermSumsAvg(residues);

  // Suffix sums for C-terminal ions
  // suffix[i] = sum of last i residues
  const monoSuffix = new Array(N + 1).fill(0);
  const avgSuffix  = new Array(N + 1).fill(0);
  for (let i = 1; i <= N; i++) {
    monoSuffix[i] = monoSuffix[i - 1] + residues[N - i].monoMass;
    avgSuffix[i]  = avgSuffix[i - 1]  + residues[N - i].avgMass;
  }

  // Determine which ion columns to render and in what order:
  // Left side (N-terminal): highest charge → lowest for each ion type
  // Right side (C-terminal): lowest charge → highest for each ion type
  // Order: a, b, c | [#] [seq] [#rev] | y, z

  const nTermTypes = [];
  if (useA) nTermTypes.push("a");
  if (useB) nTermTypes.push("b");
  if (useC) nTermTypes.push("c");

  const cTermTypes = [];
  if (useY) cTermTypes.push("y");
  if (useZ) cTermTypes.push("z");

  const charges = [];
  for (let z = fMin; z <= fMax; z++) charges.push(z);
  const chargesRev = [...charges].reverse();

  // Build header columns
  const headers = [];
  // N-terminal: for each type, highest charge first
  nTermTypes.forEach(type => {
    chargesRev.forEach(z => {
      headers.push({ side: "n", type, charge: z });
    });
  });
  // Center columns
  headers.push({ side: "center", type: "idx" });
  headers.push({ side: "center", type: "seq" });
  headers.push({ side: "center", type: "idxrev" });
  // C-terminal: for each type, lowest charge first
  cTermTypes.forEach(type => {
    charges.forEach(z => {
      headers.push({ side: "c", type, charge: z });
    });
  });

  // Header row HTML
  let html = `<table><thead><tr>`;
  headers.forEach(h => {
    if (h.side === "center") {
      const label = h.type === "idx" ? "#" : h.type === "idxrev" ? "#" : "AA";
      html += `<th class="center-col">${label}</th>`;
    } else {
      const superscript = h.charge > 1 ? `${h.charge}+` : `+`;
      html += `<th class="${h.type}-ion">${h.type}<sup>${superscript}</sup></th>`;
    }
  });
  html += `</tr></thead><tbody>`;

  // Data rows: row i corresponds to residue i (1-indexed)
  for (let i = 1; i <= N; i++) {
    html += `<tr>`;
    headers.forEach(h => {
      if (h.side === "center") {
        if (h.type === "idx")    html += `<td class="idx-cell">${i}</td>`;
        if (h.type === "seq")    html += `<td class="seq-cell">${residues[i - 1].label}</td>`;
        if (h.type === "idxrev") html += `<td class="idx-cell">${N - i + 1}</td>`;
        return;
      }

      const z = h.charge;
      const useAvg = (massType === "average" && z >= 2);

      // N-terminal ions: b1 to b(N-1); row N gets —
      if (h.side === "n") {
        if (i === N) { html += `<td class="dash">—</td>`; return; }
        const prefMono = monoPrefix[i];   // sum of first i residues
        const prefAvg  = avgPrefix[i];
        let neutral;
        if (h.type === "b") {
          neutral = useAvg ? prefAvg : prefMono;
        } else if (h.type === "a") {
          neutral = (useAvg ? prefAvg : prefMono) - CO;
        } else if (h.type === "c") {
          neutral = (useAvg ? prefAvg : prefMono) + NH3;
        }
        const mz = (neutral + z * PROTON) / z;
        html += `<td class="${h.type}-ion">${fmt(mz)}</td>`;
        return;
      }

      // C-terminal ions:
      // Row i → y ion covers residues (i+1)→N → suffix length = N - i
      // Row 1: suffLen = N-1 but this would be the near-full peptide y(N-1)... 
      //   Per requested layout: row 1 (first AA) → — for y ions
      //   Row N (last AA)  → y1 shown
      if (h.side === "c") {
        if (i === 1) { html += `<td class="dash">—</td>`; return; }  // no y for row 1
        const suffLen = N - i;   // row 2 → N-2 ... row N → 0+1 = need i-based
        // y-ion index = N - i + 1 residues from C-terminus starting AFTER position i-1
        // We want: row 2 → y(N-1), row 3 → y(N-2), ... row N → y1
        // suffix length for row i = N - i + 1... but row N must give y1 = 1 residue
        // suffLen for row i = N - (i - 1) = N - i + 1
        const sLen = N - i + 1;  // row 2→N-1 residues, row N→1 residue ✓
        const sufMono = monoSuffix[sLen];
        const sufAvg  = avgSuffix[sLen];
        let neutral;
        if (h.type === "y") {
          neutral = (useAvg ? sufAvg : sufMono) + (useAvg ? WATER_AVG : WATER_MONO);
        } else if (h.type === "z") {
          neutral = (useAvg ? sufAvg : sufMono) + (useAvg ? WATER_AVG : WATER_MONO) - NH3;
        }
        const mz = (neutral + z * PROTON) / z;
        html += `<td class="${h.type}-ion">${fmt(mz)}</td>`;
        return;
      }
    });
    html += `</tr>`;
  }

  html += `</tbody></table>`;
  return html;
}

// ─── Clear ───────────────────────────────────────────────────
clearBtn.addEventListener("click", () => {
  peptideInput.value = "";
  precursorMinEl.value = 1;
  precursorMaxEl.value = 3;
  fragmentMinEl.value  = 1;
  fragmentMaxEl.value  = 2;
  ionB.checked = true;
  ionY.checked = true;
  ionA.checked = false;
  ionC.checked = false;
  ionZ.checked = false;
  document.querySelector('input[name="mass-type"][value="mono"]').checked = true;
  document.querySelectorAll("#static-mods-container input[type=checkbox]")
    .forEach(cb => cb.checked = false);
  resultsSection.classList.add("hidden");
});

// ─── Export PDF ───────────────────────────────────────────────
exportPdfBtn.addEventListener("click", () => window.print());

// ─── Export CSV ───────────────────────────────────────────────
exportCsvBtn.addEventListener("click", exportCSV);

function exportCSV() {
  let csv = "";

  // Precursor table
  const precTable = precursorTableEl.querySelector("table");
  if (precTable) {
    csv += "PRECURSOR MASSES\n";
    csv += "Sequence," + precursorSeqEl.textContent + "\n";
    csv += tableToCSV(precTable);
    csv += "\n";
  }

  // Fragment table
  const fragTable = fragmentTableEl.querySelector("table");
  if (fragTable) {
    csv += "FRAGMENT ION MASSES\n";
    csv += tableToCSV(fragTable);
  }

  if (!csv) { alert("Nothing to export. Please calculate first."); return; }

  const blob = new Blob([csv], { type: "text/csv" });
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement("a");
  a.href     = url;
  a.download = "peptide_fragments.csv";
  a.click();
  URL.revokeObjectURL(url);
}

function tableToCSV(table) {
  let csv = "";
  table.querySelectorAll("tr").forEach(row => {
    const cells = [...row.querySelectorAll("th, td")]
      .map(c => `"${c.textContent.replace(/"/g, '""').trim()}"`);
    csv += cells.join(",") + "\n";
  });
  return csv;
}
