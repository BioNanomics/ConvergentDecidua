## Project Context

This is **ConvergentDecidua** — a Python 3.11 bioinformatics project building a comparative decidualization atlas. See [BACKGROUND.md](../BACKGROUND.md) for the scientific objective, data models, and repo structure.

### Planning — single source of truth

[PLAN.md](../PLAN.md) at the repository root is the **single source of truth** for project planning, priorities, and progress. Always read PLAN.md at the start of a planning conversation and update it (in the same PR / commit) whenever the plan changes.

**Do not** create or maintain parallel plans in:

- `/memories/session/` (workspaceStorage session memory)
- `/memories/repo/` or `/memories/` for plan content (these are for conventions, gotchas, and short-lived working notes — not the plan itself)
- Ad-hoc `PLAN_*.md`, `ROADMAP.md`, or `TODO.md` files
- Long-lived chat threads or scratchpads

If you find yourself drafting a plan in session memory, stop and mirror it into PLAN.md instead. Session/repo memory may reference PLAN.md but must not duplicate or diverge from it.

### Key conventions

- **Language**: Python 3.11, pure pip/uv (no conda). Package defined in `pyproject.toml`.
- **CLI**: `wombat` command built with Click. Entry point: `wombat/cli.py`.
- **Workflow**: Snakemake rules in `workflows/rules/`. Each rule should be callable via `wombat` CLI too.
- **Data format**: AnnData (`.h5ad`) for single-cell, Parquet for tables, DuckDB for queries.
- **Testing**: pytest. Tests in `tests/`.
- **Linting**: ruff (config in `pyproject.toml`).
- **CI**: GitHub Actions — lint, test, validate-configs, validate-workflow.

### Pre-commit verification (run before claiming "green" or committing)

Both lint **and** format must be checked. `ruff check .` alone is not
sufficient — CI runs `ruff format --check .` as a separate step and a
single unformatted file will fail the build.

```bash
ruff check .              # lint rules
ruff format --check .     # formatting — DO NOT SKIP
pytest -q                 # unit tests
pytest -m real_data -q    # real-data smoke tests (when results/ is populated)
snakemake -n --snakefile workflows/Snakefile --forceall   # DAG dry-run
```

If `ruff format --check .` reports "Would reformat", run
`ruff format <files>` and re-verify before committing. New files added
mid-session are the usual culprit — always sweep formatting on every
file touched in the change, not only on the modules you actively edited.

### Module layout

```text
wombat/          # CLI and orchestration (Click commands, config loader)
src/             # Analysis modules (ingest, metadata, qc, orthologs, cell_states, scoring, reports)
decidual_atlas/  # Streamlit visualization app
configs/         # YAML configs (datasets, species, markers)
workflows/       # Snakemake rules
```

### Code patterns

- Config loading: use `wombat.config.load_config(name)` — never hardcode paths to YAML files.
- AnnData conventions: harmonized `.obs` columns include `species`, `assay`, `cycle_stage`, `cell_type`, `cell_state`, `donor`, `sample`.
- Scoring: use `src/scoring/engine.py` generic framework. Gene sets come from `configs/markers.yaml`.
- Ortholog mapping: always go through `results/orthologs/backbone.parquet` for cross-species gene mapping.

---

## Code Quality Principles

<!-- https://github.com/mieweb/template-mieweb-opensource/blob/main/.github/copilot-instructions.md -->

### 🎯 DRY (Don't Repeat Yourself)
- **Never duplicate code**: If you find yourself copying code, extract it into a reusable function
- **Single source of truth**: Each piece of knowledge should have one authoritative representation
- **Refactor mercilessly**: When you see duplication, eliminate it immediately
- **Shared utilities**: Common patterns should be abstracted into utility functions

### 💋 KISS (Keep It Simple, Stupid)
- **Simple solutions**: Prefer the simplest solution that works
- **Avoid over-engineering**: Don't add complexity for hypothetical future needs
- **Clear naming**: Functions and variables should be self-documenting
- **Small functions**: Break down complex functions into smaller, focused ones
- **Readable code**: Code should be obvious to understand at first glance

### 🧹 Folder Philosophy
- **Clear purpose**: Every folder should have a main thing that anchors its contents.
- **No junk drawers**: Don’t leave loose files without context or explanation.
- **Explain relationships**: If it’s not elegantly obvious how files fit together, add a README or note.
- **Immediate clarity**: Opening a folder should make its organizing principle clear at a glance.

### 🔄 Refactoring Guidelines
- **Continuous improvement**: Refactor as you work, not as a separate task
- **Safe refactoring**: Always run tests before and after refactoring
- **Incremental changes**: Make small, safe changes rather than large rewrites
- **Preserve behavior**: Refactoring should not change external behavior
- **Code reviews**: All refactoring should be reviewed for correctness

### ⚰️ Dead Code Management
- **Immediate removal**: Delete unused code immediately when identified
- **Historical preservation**: Move significant dead code to `.attic/` directory with context
- **Documentation**: Include comments explaining why code was moved to attic
- **Regular cleanup**: Review and clean attic directory periodically
- **No accumulation**: Don't let dead code accumulate in active codebase

### 🌐 Testing with MCP Browser
- Use MCP browser in Playwright if available to test functionality
- **Never close the browser** after running MCP browser commands unless explicitly asked
- Let the user interact with the browser after navigation or testing
- Only use `browser_close` when the user specifically requests it

## HTML & CSS Guidelines
- **Semantic Naming**: Every `<div>` and other structural element must use a meaningful, semantic class name that clearly indicates its purpose or role within the layout.
- **CSS Simplicity**: Styles should avoid global resets or overrides that affect unrelated components or default browser behavior. Keep changes scoped and minimal.
- **SASS-First Approach**: All styles should be written in SASS (SCSS) whenever possible. Each component should have its own dedicated SASS file to promote modularity and maintainability.

## Accessibility (ARIA Labeling)

### 🎯 Interactive Elements
- **All interactive elements** (buttons, links, forms, dialogs) must include appropriate ARIA roles and labels
- **Use ARIA attributes**: Implement aria-label, aria-labelledby, and aria-describedby to provide clear, descriptive information for screen readers
- **Semantic HTML**: Use semantic HTML wherever possible to enhance accessibility

### 📢 Dynamic Content
- **Announce updates**: Ensure all dynamic content updates (modals, alerts, notifications) are announced to assistive technologies using aria-live regions
- **Maintain tab order**: Maintain logical tab order and keyboard navigation for all features
- **Visible focus**: Provide visible focus indicators for all interactive elements

## Internationalization (I18N)

### 🌍 Text and Language Support
- **Externalize text**: All user-facing text must be externalized for translation
- **Multiple languages**: Support multiple languages, including right-to-left (RTL) languages such as Arabic and Hebrew
- **Language selector**: Provide a language selector for users to choose their preferred language

### 🕐 Localization
- **Format localization**: Ensure date, time, number, and currency formats are localized based on user settings
- **UI compatibility**: Test UI layouts for text expansion and RTL compatibility
- **Unicode support**: Use Unicode throughout to support international character sets

## Documentation Preferences

### Diagrams and Visual Documentation
- **Always use Mermaid diagrams** instead of ASCII art for workflow diagrams, architecture diagrams, and flowcharts
- **Use memorable names** instead of single letters in diagrams (e.g., `Engine`, `Auth`, `Server` instead of `A`, `B`, `C`)
- Use appropriate Mermaid diagram types:
  - `graph TB` or `graph LR` for workflow architectures 
  - `flowchart TD` for process flows
  - `sequenceDiagram` for API interactions
  - `gitgraph` for branch/release strategies
- Include styling with `classDef` for better visual hierarchy
- Add descriptive comments and emojis sparingly for clarity

### Documentation Standards
- Keep documentation DRY (Don't Repeat Yourself) - reference other docs instead of duplicating
- Use clear cross-references between related documentation files
- Update the main architecture document when workflow structure changes

## Working with GitHub Actions Workflows

### Development Philosophy
- **Script-first approach**: All workflows should call scripts that can be run locally
- **Local development parity**: Developers should be able to run the exact same commands locally as CI runs
- **Simple workflows**: GitHub Actions should be thin wrappers around scripts, not contain complex logic
- **Easy debugging**: When CI fails, developers can reproduce the issue locally by running the same script

## Quick Reference

### 🪶 All Changes should be considered for Pull Request Philosophy

* **Smallest viable change**: Always make the smallest change that fully solves the problem.
* **Fewest files first**: Start with the minimal number of files required.
* **No sweeping edits**: Broad refactors or multi-module changes must be split or proposed as new components.
* **Isolated improvements**: If a change grows complex, extract it into a new function, module, or component instead of modifying multiple areas.
* **Direct requests only**: Large refactors or architectural shifts should only occur when explicitly requested.

## Command-line interface (CLI) guidelines
When running long-running or potentially hanging commands, always capture output using a unique log file per process to avoid conflicts. If the output is disposable use this pattern:
```bash
LOGFILE=$(mktemp /tmp/logfile.XXXXXX)
command 2>&1 | tee "$LOGFILE" | tail -30
```
where command is is just the main command (others could be chained with complex paramters and pipes)

This ensures real-time visibility of the last 30 lines and preserves the complete output for inspection, while remaining safe for concurrent or parallel runs.

If the output is important for debugging, ensure it is preserved and easily accessible. For example:
```bash
LOGFILE="logs/command_$(date +%Y%m%d_%H%M%S).log"
command 2>&1 | tee "$LOGFILE" | tail -30
```
This approach timestamps the log file for easy identification and prevents overwriting logs from previous runs. 

### Code Quality Checklist
- [ ] **DRY**: No code duplication - extracted reusable functions?
- [ ] **KISS**: Simplest solution that works?
- [ ] **Minimal Changes**: Smallest viable change made for PR?
- [ ] **Naming**: Self-documenting function/variable names?
- [ ] **Size**: Functions small and focused?
- [ ] **Dead Code**: Removed or archived appropriately?
- [ ] **Accessibility**: ARIA labels and semantic HTML implemented?
- [ ] **I18N**: User-facing text externalized for translation?
- [ ] **Lint**: Run linter if appropriate
- [ ] **Test**: Run tests
