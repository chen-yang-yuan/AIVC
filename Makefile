# ============================
# AIVC Makefile
# ============================

.DEFAULT_GOAL := push

.PHONY: push status log data-push data-push-all out-pull data-push-dry out-pull-dry check-branch

# Git settings
BRANCH := main
REMOTE := origin
TIMESTAMP := $(shell date "+%Y-%m-%d %H:%M:%S")

# ----------------------------
# Git automation
# ----------------------------
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
HGCC_BASE := $(HOME)/hulab/projects/AIVC

LOCAL_BASE := $(HOME)/Dropbox/AIVC
LOCAL_DATA := $(LOCAL_BASE)/data
LOCAL_OUT  := $(LOCAL_BASE)/output

REMOTE_DATA := $(HGCC_BASE)/data
REMOTE_OUT  := $(HGCC_BASE)/output

RSYNC_COMMON := -avh --progress
RSYNC_FAST := $(RSYNC_COMMON) --partial --inplace

data-push:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p '$(REMOTE_DATA)'"
	rsync $(RSYNC_FAST) \
	  "$(LOCAL_DATA)"/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/

out-pull:
	@mkdir -p "$(LOCAL_OUT)"
	rsync $(RSYNC_FAST) \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)"/ \
	  "$(LOCAL_OUT)"/

data-push-dry:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p '$(REMOTE_DATA)'"
	@set -e; \
	if ls "$(LOCAL_DATA)"/Xenium* >/dev/null 2>&1; then \
	  rsync $(RSYNC_FAST) --dry-run "$(LOCAL_DATA)"/Xenium*/ \
	    $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/ ; \
	else \
	  echo "⚠️  No Xenium* folders found under $(LOCAL_DATA). Nothing to sync."; \
	fi

out-pull-dry:
	@mkdir -p "$(LOCAL_OUT)"
	rsync $(RSYNC_FAST) --dry-run \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)"/ \
	  "$(LOCAL_OUT)"/