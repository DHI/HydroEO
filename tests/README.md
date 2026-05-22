# Testing

## Running tests

Install the development environment with all testing tools:

```sh
uv sync --all-extras
```

Then run test suites:

```sh
make test                          # all tests (disabled warnings)
uv run pytest -m unit              # fast mocked tests, no network
uv run pytest -m api_contract      # endpoint/schema checks
uv run pytest -m integration -v    # live API calls (requires credentials)
make coverage                      # HTML coverage report
```

Run a specific file or test:

```sh
uv run pytest tests/unit/test_flows.py -v
uv run pytest tests/unit/test_flows.py::test_download_pld_downloads_when_missing -v
```

## Pytest markers

Defined in [`pyproject.toml`](../pyproject.toml):

| Marker | Description |
| --- | --- |
| `unit` | Fast mocked tests; no external network dependency |
| `integration` | Live tests calling external services; require valid credentials |
| `api_contract` | Endpoint reachability and response schema checks |

## Credentials for tests

Both Earthdata variable naming schemes are needed in practice — some tests use `EDL_*`, others use `EARTHDATA_*`:

```sh
export EDL_USERNAME="your_earthdata_username"
export EDL_PASSWORD="your_earthdata_password"
export EARTHDATA_USERNAME="$EDL_USERNAME"
export EARTHDATA_PASSWORD="$EDL_PASSWORD"
export CREODIAS_USERNAME="your_creodias_username"
export CREODIAS_PASSWORD="your_creodias_password"
export HYDROWEB_API_KEY="your_hydroweb_api_key"
```

## CI/CD and GitHub setup

Automated checks run on every push and pull request via GitHub Actions (see [`.github/workflows/ci.yml`](../.github/workflows/ci.yml)):

1. **Lint** (ruff) — code style and imports
2. **Unit tests** — Ubuntu and Windows, Python 3.10 and 3.12; coverage published to Codecov
3. **Integration tests** — live API calls (only on official repo when enabled)

**For repo maintainers:** To enable live integration tests, set these GitHub secrets and variables:

*Secrets:*
- `EDL_USERNAME` — Earthdata Login username
- `EDL_PASSWORD` — Earthdata Login password
- `CREODIAS_USERNAME` — Copernicus CDSE username
- `CREODIAS_PASSWORD` — Copernicus CDSE password

*Repository variables:*
- `RUN_INTEGRATION_TESTS` = `true` — enables integration tests on push (forks and PRs skip them automatically)
