"""Tests for the wombat CLI."""

from __future__ import annotations

from click.testing import CliRunner

from wombat.cli import cli


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "Wombat" in result.output


def test_cli_subcommands_exist():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    expected = [
        "init",
        "validate-config",
        "build-registry",
        "fetch",
        "qc",
        "orthologs",
        "integrate",
        "score-decidua",
        "serve-atlas",
    ]
    for cmd in expected:
        assert cmd in result.output, f"Missing subcommand: {cmd}"
