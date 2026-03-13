const state = {
  options: [],
  sections: [],
  visibleFilter: "",
};

const elements = {
  helpSource: document.getElementById("help-source"),
  parseStatus: document.getElementById("parse-status"),
  search: document.getElementById("search"),
  inputFile: document.getElementById("input-file"),
  commandPreview: document.getElementById("command-preview"),
  sections: document.getElementById("sections"),
  summary: document.getElementById("summary"),
  reloadHelp: document.getElementById("reload-help"),
  applyHelp: document.getElementById("apply-help"),
  resetAll: document.getElementById("reset-all"),
  copyCommand: document.getElementById("copy-command"),
  optionTemplate: document.getElementById("option-template"),
};

function shellEscape(raw) {
  if (raw === "") return "''";
  if (/^[a-zA-Z0-9_./:+-]+$/.test(raw)) return raw;
  return `'${String(raw).replace(/'/g, `'"'"'`)}'`;
}

function parseChoiceLine(line) {
  const trimmed = line.trim();
  const match = trimmed.match(/^([<>=]?\s*[-\w.,+()]+)\s*:\s*(.+)$/);
  if (!match) return null;
  return { key: match[1].replace(/\s+/g, " ").trim(), label: match[2].trim() };
}

function splitDescriptionAndDefault(text) {
  const idx = text.lastIndexOf(": ");
  if (idx === -1) return { description: text.trim(), defaultValue: "" };
  return {
    description: text.slice(0, idx).trim(),
    defaultValue: text.slice(idx + 2).trim(),
  };
}

function normalizeFlag(flagToken) {
  // Convert -h(--help) to -h while preserving alt long flag in description.
  return flagToken.replace(/\(.*\)/, "");
}

function parseHelpText(helpText) {
  const lines = helpText.split(/\r?\n/);
  const options = [];
  let current = null;
  let section = "Main Options";

  for (const rawLine of lines) {
    const line = rawLine.replace(/\t/g, "    ");
    const sectionMatch = line.match(/^-----\s*(.+?)\s*-----\s*$/);
    if (sectionMatch) {
      section = sectionMatch[1].trim();
      current = null;
      continue;
    }

    const typedMatch = line.match(/^\s+(-{1,2}[A-Za-z0-9][A-Za-z0-9-]*(?:\([^)]*\))?)\s+\[([IFS])\]\s+(.+)$/);
    if (typedMatch) {
      const flagToken = typedMatch[1];
      const type = typedMatch[2];
      const rest = typedMatch[3];
      const parsed = splitDescriptionAndDefault(rest);
      current = {
        section,
        flag: normalizeFlag(flagToken),
        originalFlag: flagToken,
        type,
        description: parsed.description,
        defaultValue: parsed.defaultValue,
        details: [],
        choices: [],
        enabled: false,
        value: parsed.defaultValue,
      };
      options.push(current);
      continue;
    }

    const flagOnlyMatch = line.match(/^\s+(-{1,2}[A-Za-z0-9][A-Za-z0-9-]*(?:\([^)]*\))?)\s*:\s+(.+)$/);
    if (flagOnlyMatch) {
      const flagToken = flagOnlyMatch[1];
      current = {
        section,
        flag: normalizeFlag(flagToken),
        originalFlag: flagToken,
        type: "B",
        description: flagOnlyMatch[2].trim(),
        defaultValue: "",
        details: [],
        choices: [],
        enabled: false,
        value: "",
      };
      options.push(current);
      continue;
    }

    if (current && line.trim() !== "") {
      const choice = parseChoiceLine(line);
      if (choice) {
        current.choices.push(choice);
      }
      current.details.push(line.trim());
    }
  }

  // Keep only actual options and drop decorative duplicates.
  const unique = new Map();
  for (const opt of options) {
    if (!opt.flag.startsWith("-")) continue;
    if (!unique.has(opt.flag)) unique.set(opt.flag, opt);
  }

  return Array.from(unique.values());
}

function makeInputForOption(option) {
  const wrap = document.createElement("div");

  if (option.type === "B") {
    const span = document.createElement("span");
    span.textContent = "No argument";
    span.className = "opt-meta";
    wrap.appendChild(span);
    return wrap;
  }

  const hasSimpleChoices = option.choices.length > 1 && option.choices.every((c) => !/[<>]/.test(c.key));
  if (hasSimpleChoices) {
    const select = document.createElement("select");
    for (const c of option.choices) {
      const item = document.createElement("option");
      item.value = c.key;
      item.textContent = `${c.key} -> ${c.label}`;
      if (String(c.key) === String(option.defaultValue)) {
        item.selected = true;
      }
      select.appendChild(item);
    }
    select.addEventListener("input", () => {
      option.value = select.value;
      option.enabled = true;
      renderCommand();
    });
    wrap.appendChild(select);
    return wrap;
  }

  const input = document.createElement("input");
  input.value = option.defaultValue;
  if (option.type === "I") {
    input.type = "number";
    input.step = "1";
  } else if (option.type === "F") {
    input.type = "number";
    input.step = "any";
  } else {
    input.type = "text";
  }
  input.addEventListener("input", () => {
    option.value = input.value;
    option.enabled = true;
    renderCommand();
  });
  wrap.appendChild(input);
  return wrap;
}

function renderOptions() {
  const filter = state.visibleFilter.toLowerCase();
  elements.sections.innerHTML = "";

  const grouped = new Map();
  for (const option of state.options) {
    if (!grouped.has(option.section)) grouped.set(option.section, []);
    grouped.get(option.section).push(option);
  }

  let visibleCount = 0;
  for (const [sectionName, options] of grouped) {
    const filtered = options.filter((option) => {
      const hay = `${option.flag} ${option.description} ${option.details.join(" ")}`.toLowerCase();
      return !filter || hay.includes(filter);
    });
    if (!filtered.length) continue;
    visibleCount += filtered.length;

    const details = document.createElement("details");
    details.className = "section";
    details.open = true;

    const summary = document.createElement("summary");
    summary.textContent = `${sectionName} (${filtered.length})`;
    details.appendChild(summary);

    const body = document.createElement("div");
    body.className = "section-body";

    for (const option of filtered) {
      const fragment = elements.optionTemplate.content.cloneNode(true);
      const row = fragment.querySelector(".option-row");
      const checkbox = fragment.querySelector(".opt-enabled");
      const flag = fragment.querySelector(".flag");
      const desc = fragment.querySelector(".opt-desc");
      const meta = fragment.querySelector(".opt-meta");
      const inputWrap = fragment.querySelector(".opt-input-wrap");

      flag.textContent = option.originalFlag;
      desc.textContent = option.description;

      const extras = [];
      extras.push(`Type: ${option.type === "B" ? "flag" : option.type}`);
      if (option.defaultValue !== "") extras.push(`Default: ${option.defaultValue}`);
      if (option.details.length) extras.push(option.details.join("\n"));
      meta.textContent = extras.join("\n");

      checkbox.checked = option.enabled;
      checkbox.addEventListener("change", () => {
        option.enabled = checkbox.checked;
        renderCommand();
      });

      inputWrap.appendChild(makeInputForOption(option));

      row.dataset.flag = option.flag;
      body.appendChild(fragment);
    }

    details.appendChild(body);
    elements.sections.appendChild(details);
  }

  elements.summary.textContent = `${state.options.length} options parsed, ${visibleCount} currently shown, ${state.options.filter((o) => o.enabled).length} selected.`;
}

function renderCommand() {
  const parts = ["petar"];

  for (const opt of state.options) {
    if (!opt.enabled) continue;
    if (opt.type === "B") {
      parts.push(opt.flag);
      continue;
    }
    const value = (opt.value ?? "").toString().trim();
    if (value === "") continue;
    parts.push(opt.flag, shellEscape(value));
  }

  const inputName = elements.inputFile.value.trim();
  if (inputName) parts.push(shellEscape(inputName));

  elements.commandPreview.textContent = parts.join(" ");
  elements.summary.textContent = `${state.options.length} options parsed, ${state.options.filter((o) => o.enabled).length} selected.`;
}

function parseAndRender(helpText) {
  const parsed = parseHelpText(helpText);
  if (!parsed.length) {
    elements.parseStatus.textContent = "No options found. Paste complete `petar -h` output and parse again.";
    return;
  }

  state.options = parsed;
  state.sections = [...new Set(parsed.map((o) => o.section))];
  elements.parseStatus.textContent = `Parsed ${parsed.length} options across ${state.sections.length} sections.`;

  renderOptions();
  renderCommand();
}

async function loadHelpFile() {
  try {
    const res = await fetch("./petar-help.txt", { cache: "no-store" });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const text = await res.text();
    elements.helpSource.value = text;
    parseAndRender(text);
  } catch (err) {
    elements.parseStatus.textContent = `Failed to load petar-help.txt automatically (${err.message}). Paste text manually.`;
  }
}

function resetSelections() {
  for (const opt of state.options) {
    opt.enabled = false;
    opt.value = opt.defaultValue;
  }
  renderOptions();
  renderCommand();
}

async function copyCommand() {
  try {
    await navigator.clipboard.writeText(elements.commandPreview.textContent);
    elements.parseStatus.textContent = "Command copied to clipboard.";
  } catch (err) {
    elements.parseStatus.textContent = `Clipboard copy failed: ${err.message}`;
  }
}

elements.reloadHelp.addEventListener("click", loadHelpFile);
elements.applyHelp.addEventListener("click", () => parseAndRender(elements.helpSource.value));
elements.resetAll.addEventListener("click", resetSelections);
elements.copyCommand.addEventListener("click", copyCommand);
elements.search.addEventListener("input", () => {
  state.visibleFilter = elements.search.value;
  renderOptions();
});
elements.inputFile.addEventListener("input", renderCommand);

loadHelpFile();
