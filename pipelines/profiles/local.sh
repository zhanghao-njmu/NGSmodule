#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/profiles/local.sh — execute in the current shell.
#
# This is the default; included for completeness so the dispatch path is
# uniform. The orchestrator's fast-path skips even loading this when
# NGS_PROFILE=local, but you can still source it for testing.
###############################################################################

profile_submit_sample() {
  local pipeline="${1:?pipeline required}"
  local sample="${2:?sample required}"
  local fn="run_${pipeline}_for_sample"
  "$fn" "$sample"
}
