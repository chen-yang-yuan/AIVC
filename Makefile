# ============================
# AIVC Makefile
# ============================

.DEFAULT_GOAL := push

.PHONY: push status log check-branch \
        data-push-Xenium data-push out-pull data-push-dry out-pull-dry \
        enrich-push enrich-push-dry \
        env-check env-create env-update env-recreate

# ----------------------------
# Git automation
# ----------------------------
BRANCH := main
REMOTE := origin
TIMESTAMP := $(shell date "+%Y-%m-%d %H:%M:%S")

check-branch:
	@CURRENT=$$(git branch --show-current); \
	if [ "$$CURRENT" != "$(BRANCH)" ]; then \
	  echo "❌ On branch $$CURRENT (expected $(BRANCH))"; exit 1; \
	fi

push: check-branch
	@echo "🔍 Checking git status..."
	@git rev-parse --is-inside-work-tree >/dev/null 2>&1 || \
	  (echo "❌ Not a git repository" && exit 1)
	@echo "📌 Staging changes..."
	@git add -A
	@echo "📝 Committing changes..."
	@git diff --cached --quiet || \
	  git commit -m "Auto-commit: $(TIMESTAMP)"
	@echo "🚀 Pushing to $(REMOTE)/$(BRANCH)..."
	@git push $(REMOTE) $(BRANCH)
	@echo "✅ Done."

status:
	@git status

log:
	@git --no-pager log --oneline -5

# ============================
# Data/output sync (rsync)
# ============================

HGCC_USER := cyuan36
HGCC_HOST := hgcc.emory.edu

HGCC_BASE := ~/hulab/projects/AIVC
LOCAL_BASE := $(HOME)/Dropbox/AIVC

LOCAL_DATA := $(LOCAL_BASE)/data
LOCAL_OUT  := $(LOCAL_BASE)/output

REMOTE_DATA := $(HGCC_BASE)/data
REMOTE_OUT  := $(HGCC_BASE)/output

RSYNC_COMMON := -avh --progress
RSYNC_FAST := $(RSYNC_COMMON) --partial --inplace

# Push only Xenium* subfolders from local data -> HGCC
data-push-Xenium:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p '$(REMOTE_DATA)'"
	@set -e; \
	if ls "$(LOCAL_DATA)"/Xenium* >/dev/null 2>&1; then \
	  rsync $(RSYNC_FAST) "$(LOCAL_DATA)"/Xenium* \
	    $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/ ; \
	else \
	  echo "⚠️  No Xenium* folders found under $(LOCAL_DATA). Nothing to sync."; \
	fi

# Push ALL local data -> HGCC (use with care)
data-push:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p '$(REMOTE_DATA)'"
	rsync $(RSYNC_FAST) \
	  "$(LOCAL_DATA)"/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/

# Pull HGCC output -> local output
out-pull:
	@mkdir -p "$(LOCAL_OUT)"
	rsync $(RSYNC_FAST) \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)"/ \
	  "$(LOCAL_OUT)"/

# Push the cytoplasm-enriched gene sets (output/1_exploration/enrichment, from 1_enrichment.ipynb) -> HGCC
ENRICH_SUBDIR := 1_exploration/enrichment
enrich-push:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p '$(REMOTE_OUT)/$(ENRICH_SUBDIR)'"
	rsync $(RSYNC_FAST) \
	  "$(LOCAL_OUT)/$(ENRICH_SUBDIR)"/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)/$(ENRICH_SUBDIR)"/

enrich-push-dry:
	rsync $(RSYNC_FAST) --dry-run \
	  "$(LOCAL_OUT)/$(ENRICH_SUBDIR)"/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)/$(ENRICH_SUBDIR)"/

# Pull granule_adata.h5ad from HGCC -> local data
pull_granule_data:
	@mkdir -p "$(LOCAL_DATA)"
	rsync $(RSYNC_FAST) \
	  --include='*/' \
	  --include='processed_data/' \
	  --include='processed_data/granule_adata.h5ad' \
	  --exclude='*' \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/ \
	  "$(LOCAL_DATA)"/

# Optional: preview what would transfer (no changes)
data-push-Xenium-dry:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p '$(REMOTE_DATA)'"
	@set -e; \
	if ls "$(LOCAL_DATA)"/Xenium* >/dev/null 2>&1; then \
	  rsync $(RSYNC_FAST) --dry-run "$(LOCAL_DATA)"/Xenium* \
	    $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/ ; \
	else \
	  echo "⚠️  No Xenium* folders found under $(LOCAL_DATA). Nothing to sync."; \
	fi

data-push-dry:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p '$(REMOTE_DATA)'"
	rsync $(RSYNC_FAST) --dry-run \
	  "$(LOCAL_DATA)"/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/

out-pull-dry:
	@mkdir -p "$(LOCAL_OUT)"
	rsync $(RSYNC_FAST) --dry-run \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)"/ \
	  "$(LOCAL_OUT)"/

# ============================
# Conda environment helpers
# ============================

ENV_YAML := code/utils/env.yaml
ENV_NAME := AIVC-env

env-check:
	@conda env list | grep -q "^$(ENV_NAME)\b" && \
	  echo "✅ Found $(ENV_NAME)" || \
	  echo "❌ Conda env $(ENV_NAME) not found"

# Create env only if it does not exist
env-create:
	@conda env list | grep -q "^$(ENV_NAME)\b" && \
	  echo "ℹ️  $(ENV_NAME) already exists; skipping create." || \
	  conda env create -f $(ENV_YAML)

# Update env if it exists, otherwise create it
env-update:
	@conda env list | grep -q "^$(ENV_NAME)\b" && \
	  conda env update -f $(ENV_YAML) --prune || \
	  conda env create -f $(ENV_YAML)

# Explicitly destroy and rebuild (dangerous, but sometimes necessary)
env-recreate:
	@echo "⚠️  Recreating $(ENV_NAME) from scratch..."
	conda remove -n $(ENV_NAME) --all -y
	conda env create -f $(ENV_YAML)