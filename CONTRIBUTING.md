# Contributing to HydroEO

Thanks for your interest in improving HydroEO! Bug reports, feature requests, and pull requests are welcome.

## Reporting issues

Open a [GitHub issue](https://github.com/DHI/HydroEO/issues) and include:

- what you were trying to do (scenario/config, with credentials removed),
- what happened vs. what you expected, including the full traceback or relevant log lines (`logs/HydroEO_*.log`),
- your OS, Python version, and HydroEO version.

**Never paste credentials, API keys, or tokens** into issues, configs, notebooks, or commits. Use environment variables instead (see the README).

## Development setup

Requires Python 3.11+, [uv](https://docs.astral.sh/uv/), and git.

```bash
git clone https://github.com/DHI/HydroEO.git
cd HydroEO
uv sync --all-extras
```

## Making changes

1. Create a branch from `main` (e.g. `fix/pld-download-timeout`, `feature/river-cleaning`).
2. Make your change, with tests.
3. Run the full check suite locally before pushing:

   ```bash
   make check   # ruff lint + mypy + pytest + doctests
   make format  # auto-format with ruff
   ```

4. Open a pull request against `main` describing what changed and why. Link any related issue.

CI runs ruff and the unit tests on Ubuntu and Windows (Python 3.11 and 3.13); it must pass before merging. Code owners (see `.github/CODEOWNERS`) are requested for review automatically.

## Coding standards

- **Style:** formatted and linted with ruff. Type hints use modern syntax (`str | None`).
- **Docstrings:** NumPy style for public functions and classes.
- **Logging:** use a module-level `logger = logging.getLogger(__name__)`. Don't use `print`.
- **Validation:** config validation in `Project.initialize()` collects *all* errors before doing any I/O. Keep it that way when you add new checks.
- **Configs:** if you change config keys or output structure, update both `configs/*.yaml` and the matching `configs/*.md`.

## Tests

Mark every test with one of the markers declared in `pyproject.toml` (unknown markers fail under `--strict-markers`):

| Marker | Meaning |
|---|---|
| `unit` | fast, mocked, no network — required for all new logic |
| `integration` | hits live APIs, needs credentials |
| `api_contract` | credential-free checks against external API contracts |
| `sliderule` | depends on the SlideRule API |

```bash
uv run pytest -m unit                   # what CI runs
uv run pytest -m integration -v         # needs credentials, see README
```

Shared fixtures (a mock `Project` built with `SimpleNamespace`, GeoDataFrame builders) are in `tests/conftest.py` and `tests/unit/test_flows.py`.

## License

By contributing, you agree that your contributions will be licensed under the project's [LICENSE](LICENSE).
