# ============================
# AIVC Makefile
# ============================

.DEFAULT_GOAL := push

# Git settings
BRANCH := main
REMOTE := origin

# Timestamp for auto-commit message
TIMESTAMP := $(shell date "+%Y-%m-%d %H:%M:%S")

# ----------------------------
# Main target
# ----------------------------
push:
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

# ----------------------------
# Optional helpers
# ----------------------------
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

LOCAL_BASE := ~/Dropbox/AIVC
LOCAL_DATA := $(LOCAL_BASE)/data
LOCAL_OUT  := $(LOCAL_BASE)/output

REMOTE_DATA := $(HGCC_BASE)/data
REMOTE_OUT  := $(HGCC_BASE)/output

RSYNC_COMMON := -avh --progress
# Good defaults for big folders:
RSYNC_FAST := $(RSYNC_COMMON) --partial --inplace

# Push only Xenium* subfolders from local data -> HGCC
data-push:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p $(REMOTE_DATA)"
	rsync $(RSYNC_FAST) \
	  "$(LOCAL_DATA)"/Xenium*/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/

# Push ALL local data -> HGCC (use with care)
data-push-all:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p $(REMOTE_DATA)"
	rsync $(RSYNC_FAST) \
	  "$(LOCAL_DATA)"/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/

# Pull HGCC output -> local output
out-pull:
	@mkdir -p "$(LOCAL_OUT)"
	rsync $(RSYNC_FAST) \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)"/ \
	  "$(LOCAL_OUT)"/

# Optional: preview what would transfer (no changes)
data-push-dry:
	@ssh $(HGCC_USER)@$(HGCC_HOST) "mkdir -p $(REMOTE_DATA)"
	rsync $(RSYNC_FAST) --dry-run \
	  "$(LOCAL_DATA)"/Xenium*/ \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_DATA)"/

out-pull-dry:
	@mkdir -p "$(LOCAL_OUT)"
	rsync $(RSYNC_FAST) --dry-run \
	  $(HGCC_USER)@$(HGCC_HOST):"$(REMOTE_OUT)"/ \
	  "$(LOCAL_OUT)"/