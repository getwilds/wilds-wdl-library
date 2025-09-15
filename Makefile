lint_sprocket:
	@echo "Running sprocket lint on all modules..."
	@echo "Checking if sprocket is available..."
	@if ! command -v sprocket >/dev/null 2>&1; then \
		echo >&2 "Error: sprocket is not installed or not in PATH. Install sprocket."; \
		exit 1; \
	fi
	@for dir in modules/*/; do \
		if [ -d "$$dir" ]; then \
			echo "Linting $$dir"; \
			sprocket lint "$$dir"; \
		fi; \
	done

lint_miniwdl:
	@echo "Running miniwdl lint on all modules..."
	@echo "Checking if uv is available..."
	@if ! command -v uv >/dev/null 2>&1; then \
		echo >&2 "Error: uv is not installed or not in PATH. Install uv."; \
		exit 1; \
	fi
	@for file in modules/*/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Linting $$file"; \
			uv run --python 3.13.5 --with miniwdl miniwdl check "$$file"; \
		fi; \
	done

lint: lint_sprocket lint_miniwdl
