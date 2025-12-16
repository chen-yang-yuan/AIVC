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