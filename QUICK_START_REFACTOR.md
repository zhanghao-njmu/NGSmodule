# NGSmodule 快速重构指南

## 🎯 目标
在1-2周内快速提升代码质量，修复关键bug，为Web平台开发打好基础。

---

## 📦 第一周：紧急修复和基础改进

### Day 1: 全局错误处理

#### 1.1 创建错误处理库
```bash
# 创建文件: lib/error_handler.sh
mkdir -p lib

cat > lib/error_handler.sh << 'EOF'
#!/usr/bin/env bash

# 启用严格模式
set -euo pipefail

# 错误处理器
error_handler() {
    local exit_code=$?
    local line_number=$1
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2
    echo "❌ ERROR: Command failed" >&2
    echo "   Exit Code: $exit_code" >&2
    echo "   Line: $line_number" >&2
    echo "   Script: ${BASH_SOURCE[1]}" >&2
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2

    # 打印调用栈
    echo "Stack trace:" >&2
    local frame=0
    while caller $frame >&2; do
        ((frame++)) || true
    done

    exit $exit_code
}

# 清理函数
cleanup() {
    echo "🧹 Cleaning up..." >&2
    # 清理临时文件
    [[ -n "${TMPDIR:-}" ]] && rm -rf "$TMPDIR"
    # 关闭文件描述符
    [[ -e /proc/self/fd/1000 ]] && exec 1000>&- || true
}

# 设置陷阱
trap 'error_handler $LINENO' ERR
trap 'cleanup; exit 130' INT TERM
trap 'cleanup' EXIT

# 导出函数
export -f error_handler
export -f cleanup
EOF

chmod +x lib/error_handler.sh
```

#### 1.2 在主脚本中应用
```bash
# 在 NGSmodule 文件顶部添加（第一行之后）
source "$(dirname "$0")/lib/error_handler.sh"
```

#### 1.3 批量更新所有脚本
```bash
# 创建脚本: scripts/add_error_handler.sh
cat > scripts/add_error_handler.sh << 'EOF'
#!/usr/bin/env bash

# 为所有Shell脚本添加错误处理
find . -name "*.sh" -type f | while read -r script; do
    # 跳过已经有的
    if grep -q "lib/error_handler.sh" "$script"; then
        echo "Skip: $script (already has error handler)"
        continue
    fi

    # 备份
    cp "$script" "$script.bak"

    # 在shebang后添加source语句
    sed -i '2i\\n# Error handling\nsource "$(dirname "$0")/lib/error_handler.sh"\n' "$script"

    echo "✓ Updated: $script"
done
EOF

chmod +x scripts/add_error_handler.sh
./scripts/add_error_handler.sh
```

---

### Day 2: 参数化配置

#### 2.1 创建配置模板
```bash
# 创建文件: config/defaults.sh
mkdir -p config

cat > config/defaults.sh << 'EOF'
#!/usr/bin/env bash

# NGSmodule 默认配置
# 可通过环境变量覆盖

# 参考数据路径
: "${IGENOMES_DIR:=/data/reference/iGenomes}"
: "${FASTQ_SCREEN_CONFIG:=/data/reference/FastQ_Screen/FastQ_Screen_Genomes/fastq_screen.conf}"
: "${SORTMERNA_DIR:=/data/reference/SortmeRNA}"

# 系统资源
: "${THREADS:=$(nproc)}"
: "${MAX_MEMORY:=32G}"

# 临时目录
: "${TMPDIR:=/tmp/ngsmodule.$$}"
mkdir -p "$TMPDIR"

# 日志设置
: "${LOG_LEVEL:=INFO}"  # DEBUG, INFO, WARNING, ERROR
: "${LOG_DIR:=$PWD/logs}"
mkdir -p "$LOG_DIR"

# 导出所有变量
export IGENOMES_DIR FASTQ_SCREEN_CONFIG SORTMERNA_DIR
export THREADS MAX_MEMORY TMPDIR LOG_LEVEL LOG_DIR

# 打印配置（仅在DEBUG模式）
if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
    echo "=== NGSmodule Configuration ===" >&2
    echo "IGENOMES_DIR: $IGENOMES_DIR" >&2
    echo "THREADS: $THREADS" >&2
    echo "TMPDIR: $TMPDIR" >&2
    echo "LOG_DIR: $LOG_DIR" >&2
    echo "===============================" >&2
fi
EOF
```

#### 2.2 创建用户配置示例
```bash
# 创建文件: config/user.config.example
cat > config/user.config.example << 'EOF'
# NGSmodule 用户配置示例
# 复制此文件为 ~/.ngsmodule/config.sh 并修改

# 参考数据路径（根据你的实际路径修改）
export IGENOMES_DIR="/path/to/your/iGenomes"
export FASTQ_SCREEN_CONFIG="/path/to/your/fastq_screen.conf"
export SORTMERNA_DIR="/path/to/your/SortmeRNA"

# 系统资源（根据你的服务器配置）
export THREADS=32
export MAX_MEMORY="64G"

# 临时目录（确保有足够空间）
export TMPDIR="/scratch/ngsmodule"

# 日志级别
export LOG_LEVEL="INFO"  # DEBUG, INFO, WARNING, ERROR
EOF
```

#### 2.3 更新LoadConfig.sh
```bash
# 在 LoadConfig.sh 开头添加
# 加载默认配置
source "$(dirname "$0")/config/defaults.sh"

# 加载用户配置（如果存在）
[[ -f "$HOME/.ngsmodule/config.sh" ]] && source "$HOME/.ngsmodule/config.sh"

# 使用配置变量替换硬编码路径
# 例如：将 /data/reference/iGenomes 替换为 $IGENOMES_DIR
```

---

### Day 3: 统一日志系统

#### 3.1 创建日志库
```bash
# 创建文件: lib/logging.sh
cat > lib/logging.sh << 'EOF'
#!/usr/bin/env bash

# 颜色定义
readonly COLOR_RESET="\033[0m"
readonly COLOR_RED="\033[0;31m"
readonly COLOR_GREEN="\033[0;32m"
readonly COLOR_YELLOW="\033[0;33m"
readonly COLOR_BLUE="\033[0;34m"
readonly COLOR_MAGENTA="\033[0;35m"
readonly COLOR_CYAN="\033[0;36m"

# 日志级别
declare -A LOG_LEVELS=([DEBUG]=0 [INFO]=1 [WARNING]=2 [ERROR]=3)
: "${LOG_LEVEL:=INFO}"

# 日志文件
: "${LOG_FILE:=/dev/stderr}"

# 通用日志函数
_log() {
    local level=$1
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    # 检查日志级别
    if (( ${LOG_LEVELS[$level]:-1} >= ${LOG_LEVELS[$LOG_LEVEL]:-1} )); then
        # 选择颜色
        local color=$COLOR_RESET
        case $level in
            ERROR)   color=$COLOR_RED ;;
            WARNING) color=$COLOR_YELLOW ;;
            INFO)    color=$COLOR_GREEN ;;
            DEBUG)   color=$COLOR_CYAN ;;
        esac

        # 输出日志
        echo -e "${color}[$timestamp] [$level]${COLOR_RESET} $message" >&2

        # 同时写入日志文件
        if [[ "$LOG_FILE" != "/dev/stderr" ]]; then
            echo "[$timestamp] [$level] $message" >> "$LOG_FILE"
        fi
    fi
}

# 便捷函数
log_debug()   { _log DEBUG "$@"; }
log_info()    { _log INFO "$@"; }
log_warning() { _log WARNING "$@"; }
log_error()   { _log ERROR "$@"; }
log_success() { echo -e "${COLOR_GREEN}✓${COLOR_RESET} $*" >&2; }
log_fail()    { echo -e "${COLOR_RED}✗${COLOR_RESET} $*" >&2; }

# 执行命令并记录
log_command() {
    local cmd="$*"
    log_debug "Executing: $cmd"

    local output
    local exit_code

    # 执行命令并捕获输出
    if output=$(eval "$cmd" 2>&1); then
        exit_code=0
        log_debug "Command succeeded"
        [[ -n "$output" ]] && log_debug "Output: $output"
    else
        exit_code=$?
        log_error "Command failed with exit code $exit_code"
        log_error "Command: $cmd"
        [[ -n "$output" ]] && log_error "Output: $output"
    fi

    return $exit_code
}

# 进度条
log_progress() {
    local current=$1
    local total=$2
    local label=${3:-"Progress"}

    local percent=$((current * 100 / total))
    local bar_length=50
    local filled=$((current * bar_length / total))
    local empty=$((bar_length - filled))

    printf "\r${COLOR_CYAN}$label${COLOR_RESET} ["
    printf "%${filled}s" | tr ' ' '='
    printf "%${empty}s" | tr ' ' '-'
    printf "] $percent%% ($current/$total)"

    if [[ $current -eq $total ]]; then
        printf "\n"
    fi
}

# 导出函数
export -f _log log_debug log_info log_warning log_error
export -f log_success log_fail log_command log_progress
EOF

chmod +x lib/logging.sh
```

#### 3.2 替换现有日志调用
```bash
# 在脚本中
source "$(dirname "$0")/lib/logging.sh"

# 替换:
# echo "INFO: ..." → log_info "..."
# echo "ERROR: ..." → log_error "..."
# echo "WARNING: ..." → log_warning "..."
```

---

### Day 4-5: 修复Shell脚本问题

#### 4.1 添加管道错误检查
```bash
# 创建脚本: scripts/fix_pipelines.sh
cat > scripts/fix_pipelines.sh << 'EOF'
#!/usr/bin/env bash

# 在所有管道前添加 set -o pipefail
find . -name "*.sh" -type f | while read -r script; do
    # 检查是否有管道
    if grep -q "|" "$script"; then
        # 备份
        cp "$script" "$script.bak"

        # 在管道前添加 set -o pipefail
        # 这里需要手动检查和修改
        echo "Check pipelines in: $script"
    fi
done
EOF
```

#### 4.2 修复变量引用
```bash
# 创建脚本: scripts/fix_quotes.sh
cat > scripts/fix_quotes.sh << 'EOF'
#!/usr/bin/env bash

# 常见的未引用变量模式
find . -name "*.sh" -type f | while read -r script; do
    echo "Checking: $script"

    # 检查常见问题
    grep -n 'rm -f $[A-Za-z_]' "$script" && \
        echo "  ⚠️  Warning: unquoted variable in rm command"

    grep -n 'cd $[A-Za-z_]' "$script" && \
        echo "  ⚠️  Warning: unquoted variable in cd command"
done
EOF

chmod +x scripts/fix_quotes.sh
./scripts/fix_quotes.sh
```

---

### Day 6-7: R脚本改进

#### 6.1 创建R包检查工具
```r
# 创建文件: lib/check_packages.R

check_and_install_packages <- function(required_packages) {
  #' Check and Install Required R Packages
  #'
  #' @param required_packages Named list of package:version pairs
  #' @return TRUE if all packages are available, stops otherwise

  missing_packages <- character()
  outdated_packages <- character()

  for (pkg in names(required_packages)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    } else {
      installed_version <- as.character(packageVersion(pkg))
      required_version <- required_packages[[pkg]]

      if (compareVersion(installed_version, required_version) < 0) {
        outdated_packages <- c(outdated_packages,
                               sprintf("%s (installed: %s, required: %s)",
                                       pkg, installed_version, required_version))
      }
    }
  }

  # 报告结果
  if (length(missing_packages) > 0) {
    stop(sprintf(
      "Missing required packages:\n  %s\n\nPlease install using:\n  BiocManager::install(c(%s))",
      paste(missing_packages, collapse="\n  "),
      paste(sprintf('"%s"', missing_packages), collapse=", ")
    ))
  }

  if (length(outdated_packages) > 0) {
    warning(sprintf(
      "Outdated packages detected:\n  %s\n\nConsider updating them.",
      paste(outdated_packages, collapse="\n  ")
    ))
  }

  message("✓ All required packages are available")
  invisible(TRUE)
}

# 定义所有依赖
ngs_dependencies <- list(
  # 数据处理
  "dplyr" = "1.0.0",
  "data.table" = "1.14.0",
  "tidyr" = "1.1.0",

  # 生物信息学
  "limma" = "3.50.0",
  "edgeR" = "3.36.0",
  "DESeq2" = "1.34.0",
  "Rsubread" = "2.8.0",

  # 可视化
  "ggplot2" = "3.3.0",
  "ComplexHeatmap" = "2.10.0",

  # 单细胞
  "Seurat" = "4.3.0",
  "SingleCellExperiment" = "1.16.0"
)

# 在脚本开头调用
if (interactive()) {
  check_and_install_packages(ngs_dependencies)
}
```

#### 6.2 统一R错误处理
```r
# 创建文件: lib/error_handler.R

#' Safe File Reading with Validation
safe_read_csv <- function(file_path, required_cols = NULL) {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }

  if (file.size(file_path) == 0) {
    stop(sprintf("File is empty: %s", file_path))
  }

  tryCatch({
    df <- read.csv(file_path, stringsAsFactors = FALSE, header = TRUE)

    if (nrow(df) == 0) {
      stop(sprintf("File contains no data: %s", file_path))
    }

    # 检查必需列
    if (!is.null(required_cols)) {
      missing_cols <- setdiff(required_cols, colnames(df))
      if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns: %s",
                     paste(missing_cols, collapse = ", ")))
      }
    }

    return(df)
  }, error = function(e) {
    stop(sprintf("Error reading file %s: %s", file_path, e$message))
  })
}

#' Safe File Writing with Atomic Operations
safe_write_table <- function(data, file_path, ...) {
  # 写入临时文件
  temp_file <- paste0(file_path, ".tmp")

  tryCatch({
    write.table(data, file = temp_file, ...)

    # 原子性重命名
    file.rename(temp_file, file_path)

    message(sprintf("✓ Successfully wrote: %s", file_path))
  }, error = function(e) {
    # 清理临时文件
    if (file.exists(temp_file)) {
      file.remove(temp_file)
    }
    stop(sprintf("Error writing file %s: %s", file_path, e$message))
  })
}

#' Execute with Progress Reporting
with_progress <- function(expr, message = "Processing") {
  cat(sprintf("%s...\n", message))
  start_time <- Sys.time()

  result <- tryCatch({
    expr
  }, error = function(e) {
    cat(sprintf("✗ Failed: %s\n", e$message))
    stop(e)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("✓ Completed in %.2f seconds\n", elapsed))

  invisible(result)
}
```

---

## 📦 第二周：功能增强和测试

### Day 8-9: 断点续传功能

```bash
# 创建文件: lib/checkpoint.sh
cat > lib/checkpoint.sh << 'EOF'
#!/usr/bin/env bash

: "${CHECKPOINT_DIR:=.checkpoints}"
mkdir -p "$CHECKPOINT_DIR"

# 保存检查点
checkpoint_save() {
    local sample=$1
    local step=$2
    local status=${3:-completed}
    local timestamp=$(date +%s)

    local checkpoint_file="$CHECKPOINT_DIR/${sample}.checkpoint"
    echo "${step}|${status}|${timestamp}" >> "$checkpoint_file"

    log_debug "Checkpoint saved: $sample - $step - $status"
}

# 加载检查点
checkpoint_load() {
    local sample=$1
    local step=$2
    local checkpoint_file="$CHECKPOINT_DIR/${sample}.checkpoint"

    if [[ -f "$checkpoint_file" ]]; then
        grep "^${step}|completed|" "$checkpoint_file" &>/dev/null
        return $?
    fi

    return 1
}

# 检查是否应该跳过
checkpoint_should_skip() {
    local sample=$1
    local step=$2
    local force=${3:-false}

    if [[ "$force" == "true" ]]; then
        log_info "Force mode: not skipping $step for $sample"
        return 1  # 不跳过
    fi

    if checkpoint_load "$sample" "$step"; then
        log_info "Skipping $step for $sample (already completed)"
        return 0  # 跳过
    fi

    return 1  # 不跳过
}

# 清除检查点
checkpoint_clear() {
    local sample=${1:-""}

    if [[ -n "$sample" ]]; then
        rm -f "$CHECKPOINT_DIR/${sample}.checkpoint"
        log_info "Cleared checkpoints for: $sample"
    else
        rm -rf "$CHECKPOINT_DIR"
        mkdir -p "$CHECKPOINT_DIR"
        log_info "Cleared all checkpoints"
    fi
}

export -f checkpoint_save checkpoint_load checkpoint_should_skip checkpoint_clear
EOF
```

### Day 10: 资源监控

```bash
# 创建文件: lib/resource_monitor.sh
cat > lib/resource_monitor.sh << 'EOF'
#!/usr/bin/env bash

# 启动资源监控
monitor_start() {
    local interval=${1:-60}
    local log_file=${2:-resource_usage.log}

    # 后台运行监控
    {
        echo "timestamp,cpu_percent,mem_percent,disk_percent" > "$log_file"

        while true; do
            local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
            local cpu=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)
            local mem=$(free | grep Mem | awk '{printf("%.2f", $3/$2 * 100.0)}')
            local disk=$(df -h "$PWD" | tail -1 | awk '{print $5}' | cut -d'%' -f1)

            echo "$timestamp,$cpu,$mem,$disk" >> "$log_file"
            sleep "$interval"
        done
    } &

    echo $! > "${log_file}.pid"
    log_info "Resource monitoring started (PID: $(cat ${log_file}.pid))"
}

# 停止资源监控
monitor_stop() {
    local log_file=${1:-resource_usage.log}
    local pid_file="${log_file}.pid"

    if [[ -f "$pid_file" ]]; then
        local pid=$(cat "$pid_file")
        kill "$pid" 2>/dev/null || true
        rm -f "$pid_file"
        log_info "Resource monitoring stopped"
    fi
}

export -f monitor_start monitor_stop
EOF
```

### Day 11-12: 测试框架

```bash
# 创建目录结构
mkdir -p tests/{unit,integration,fixtures}

# 创建测试框架
cat > tests/test_framework.sh << 'EOF'
#!/usr/bin/env bash

# 测试统计
TEST_COUNT=0
TEST_PASSED=0
TEST_FAILED=0
FAILED_TESTS=()

# 断言函数
assert_equals() {
    local expected=$1
    local actual=$2
    local message=${3:-""}

    ((TEST_COUNT++))

    if [[ "$expected" == "$actual" ]]; then
        ((TEST_PASSED++))
        echo "  ✓ Test $TEST_COUNT: $message"
        return 0
    else
        ((TEST_FAILED++))
        echo "  ✗ Test $TEST_COUNT: $message"
        echo "    Expected: $expected"
        echo "    Actual: $actual"
        FAILED_TESTS+=("$message")
        return 1
    fi
}

assert_file_exists() {
    local file=$1
    local message=${2:-"File $file should exist"}

    ((TEST_COUNT++))

    if [[ -f "$file" ]]; then
        ((TEST_PASSED++))
        echo "  ✓ Test $TEST_COUNT: $message"
        return 0
    else
        ((TEST_FAILED++))
        echo "  ✗ Test $TEST_COUNT: $message"
        FAILED_TESTS+=("$message")
        return 1
    fi
}

assert_command_succeeds() {
    local cmd="$*"
    local message="Command should succeed: $cmd"

    ((TEST_COUNT++))

    if eval "$cmd" &>/dev/null; then
        ((TEST_PASSED++))
        echo "  ✓ Test $TEST_COUNT: $message"
        return 0
    else
        ((TEST_FAILED++))
        echo "  ✗ Test $TEST_COUNT: $message"
        FAILED_TESTS+=("$message")
        return 1
    fi
}

# 运行测试套件
run_test_suite() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Running Test Suite"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # 运行所有测试文件
    for test_file in tests/unit/test_*.sh; do
        if [[ -f "$test_file" ]]; then
            echo ""
            echo "Running: $(basename $test_file)"
            source "$test_file"
        fi
    done

    # 总结
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Test Results"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Total: $TEST_COUNT"
    echo "Passed: $TEST_PASSED"
    echo "Failed: $TEST_FAILED"

    if (( TEST_FAILED > 0 )); then
        echo ""
        echo "Failed tests:"
        for test in "${FAILED_TESTS[@]}"; do
            echo "  - $test"
        done
        exit 1
    else
        echo ""
        echo "✓ All tests passed!"
        exit 0
    fi
}

export -f assert_equals assert_file_exists assert_command_succeeds run_test_suite
EOF

# 创建示例测试
cat > tests/unit/test_logging.sh << 'EOF'
#!/usr/bin/env bash

source "$(dirname "$0")/../test_framework.sh"
source "$(dirname "$0")/../../lib/logging.sh"

# 测试日志函数
test_log_functions() {
    # 测试log_info
    output=$(log_info "test message" 2>&1)
    assert_equals 0 $? "log_info should succeed"

    # 测试log_error
    output=$(log_error "error message" 2>&1)
    assert_equals 0 $? "log_error should succeed"
}

test_log_functions
EOF

chmod +x tests/test_framework.sh
chmod +x tests/unit/test_logging.sh
```

### Day 13-14: 文档编写

```bash
# 创建文档目录
mkdir -p docs/{user-guide,developer-guide,api}

# 创建安装文档
cat > docs/installation.md << 'EOF'
# NGSmodule 安装指南

## 系统要求

- Linux操作系统 (Ubuntu 20.04+ / CentOS 7+)
- Bash 4.0+
- Conda/Miniconda
- 至少100GB磁盘空间
- 建议16GB+内存

## 快速安装

### 1. 克隆仓库
```bash
git clone https://github.com/your-org/NGSmodule.git
cd NGSmodule
```

### 2. 配置环境
```bash
# 复制配置模板
mkdir -p ~/.ngsmodule
cp config/user.config.example ~/.ngsmodule/config.sh

# 编辑配置文件
vim ~/.ngsmodule/config.sh
```

### 3. 检查环境
```bash
./CheckENV.sh
```

### 4. 运行测试
```bash
./tests/test_framework.sh
```

## 详细配置

参见 [配置指南](user-guide/configuration.md)
EOF

# 创建快速开始文档
cat > docs/quickstart.md << 'EOF'
# 快速开始

## RNA-seq分析示例

### 1. 准备配置文件
```bash
./PreparationSteps/CreateConfigFile.sh
```

### 2. 准备样本信息
```bash
./PreparationSteps/CreateSampleInfoFile.sh
```

### 3. 创建工作目录
```bash
./PreparationSteps/CreateWorkDir.sh \\
  --maindir /path/to/project \\
  --SampleInfoFile sample_info.csv
```

### 4. 运行分析
```bash
./NGSmodule \\
  --maindir /path/to/project \\
  --step preAlignmentQC \\
  --force FALSE
```

## 更多示例

参见 [用户指南](user-guide/)
EOF
```

---

## ✅ 完成检查清单

### 第一周
- [ ] 创建lib/目录和基础库文件
- [ ] 创建config/目录和配置模板
- [ ] 创建scripts/目录和辅助脚本
- [ ] 添加全局错误处理到所有脚本
- [ ] 参数化所有硬编码路径
- [ ] 统一日志系统
- [ ] 修复Shell脚本变量引用问题
- [ ] 改进R脚本错误处理

### 第二周
- [ ] 实现断点续传功能
- [ ] 添加资源监控
- [ ] 创建测试框架
- [ ] 编写单元测试
- [ ] 创建文档结构
- [ ] 编写基础文档
- [ ] 运行全面测试
- [ ] 修复发现的问题

---

## 🎯 验证成功标准

运行以下命令确认改进成功：

```bash
# 1. 测试套件通过
./tests/test_framework.sh

# 2. ShellCheck检查
shellcheck -S warning **/*.sh

# 3. 运行小规模测试数据
./NGSmodule --maindir test_data --step preAlignmentQC

# 4. 检查日志
cat logs/*.log

# 5. 验证配置
source config/defaults.sh
env | grep -E "IGENOMES|THREADS"
```

---

## 📞 获取帮助

- 详细计划: `DEVELOPMENT_PLAN.md`
- 问题汇总: 见Task分析报告
- 架构文档: `DEVELOPMENT_PLAN.md` 第2-8节

---

**祝重构顺利！** 🚀
