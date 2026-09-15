# KFV / FGO Python

工作目录：`KFV-FGO-Comparison-github/kfv_fgo/`。基于上游 `python_colab` 分支的 `1310e8d`；外层原版 core、data、config 和 MATLAB 文件保持不变，旧入口不依赖本子目录。

## 代码结构

```text
kfv_fgo/
├── core/
│   ├── estimator.py       # KfvEstimator、FgoEstimator
│   ├── filter.py          # EKF、iEKF、rEKF、riEKF
│   ├── fgo.py             # 状态、因子、因子图、求解、边缘化
│   └── evaluation.py      # 真值评估、绘图、结果保存
├── model/
│   ├── simulation.py      # 四维仿真传播和观测
│   └── gnss_ins.py        # 十维 GNSS/INS 模型与坐标转换
├── data/                  # MAT/RINEX/IMU 读取、RTKLIB 适配、默认数据
├── config/                # settings.py 和四份 JSON
├── examples/              # 四个实验＋一个结果比较脚本
├── tests/
├── results/               # 运行生成
├── Colab.ipynb
├── pyproject.toml
└── README.md
```

core 中每类职责只有一个实现文件，不再按类建立子文件夹。示例直接读取配置和数据、调用估计器、评估保存。使用正常 Python 导入，不使用动态模块别名、转发文件或 bootstrap。

```python
from kfv_fgo.core.estimator import KfvEstimator, FgoEstimator
from kfv_fgo.config.settings import load_config
```

FGO 模板保留原分支：`imitate_kfv=True` 且窗口为 1 时直接返回 KFV。SWFGO 设置 `imitate_kfv=False`，独立执行因子图，保存完整历史。模板与 KFV 相等反映调用关系，不是两个独立 Python 求解器的证明。

## 安装

在本目录安装。验证环境为 Python 3.12，真实 GNSS/INS 使用 `pyrtklib==0.2.7`：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[real,test]"
```

下文的 `python` 指该环境的解释器。仅仿真或结果文件比较可安装 `.[plots,test]`，不需要 pyRTKLIB；关闭绘图时基本依赖 `.` 即可。安装后也可从其他当前目录通过绝对路径运行示例。

## 四个实验

```powershell
python examples/example_sw_fgo_simulation.py --window 5
python examples/example_sw_fgo_real.py --window 5
python examples/example_kfv_fgo_simulation.py
python examples/example_kfv_fgo_real.py
```

SWFGO 仿真默认窗口 1、迭代 1 次、无鲁棒核；实测默认窗口 5、迭代 5 次、无鲁棒核。可用 `--window`、`--iterations`、`--kernel` 调整。KFV/FGO 对比默认运行四种模式，可用 `--modes EKF riEKF` 选择。

所有实验支持 `--config`、`--output`、`--no-plots`，实测支持 `--data-dir`。默认数据在 data 内，无需其他代码库；urban_nav_deep 有 438 个历元，真值评估在时间重叠范围内覆盖 437 个历元，不外推。

每个示例提供可直接调用的 `run(...)` 函数，返回包含完整 X、timestamps、debug_info、配置和耗时的结果，同时保存 CSV、NPZ、指标 JSON 和图。默认输出位于 results 下与示例对应的目录；保存多次实验时指定不同的 --output。

## MATLAB/Python 结果文件比较

仿真与实测共用 **examples/compare_matlab_results.py**。它只读取已保存的 MATLAB/Python 结果，不启动 MATLAB，也不运行任何估计器；普通 Python 和 Colab 都能使用。

```powershell
python examples/example_kfv_fgo_simulation.py
python examples/compare_matlab_results.py --matlab tests/fixtures/matlab_simulation.mat --matlab-key reference.runs.KFV_EKF --python results/kfv_fgo_simulation/KFV_EKF.npz --output results/compare_simulation
```

真实数据使用 `tests/fixtures/matlab_real.mat`，先运行 `examples/example_kfv_fgo_real.py`，再将比较命令的 Python 文件改为 `results/kfv_fgo_real/KFV_EKF.npz`。参考文件覆盖的算法、参数和来源见 `tests/fixtures/README.md`；真实 SWFGO 参考使用 Huber 核，复现时应显式传入 `--window 5 --iterations 5 --kernel huber`。自己的 MATLAB 文件仍可通过 `--matlab` 指定。

支持两种已有 MATLAB 结构：原版 `result_ekf`、`result_fg_ekf` 等包含 X 的结果；此前参考导出的 `reference.runs.KFV_EKF` 等嵌套结果。通过 --matlab-key 选择；一个文件有多组结果时，错误提示会列出可选字段。Python 输入支持本库保存的 NPZ，也支持 MAT。

- X 使用 状态维数×历元，仿真 4 维、实测 10 维，不自动转置或截断。
- 优先读取 timestamps。仿真缺少时间时从 data.dt/settings.dt 或显式 --dt 构造，默认起点 0，可用 --start-time 指定。
- 实测必须有同一时间尺度的数值时间戳；若在其他字段，用 --matlab-time-key 指定，例如 epoch_times。不给实测虚构等间隔时间。
- 支持标准 MAT 文件；MAT v7.3 尚不支持，可将数值结果按 -v7 保存。
- 逐分量默认绝对容差：仿真 1e-8、实测 1e-6，--atol 可调整。输出不一致时仍保存报告并返回退出码 2；输入错误返回 1。

输出 comparison.csv、comparison.npz、metrics.json、轨迹图和逐状态差异图，记录输入哈希、时间来源、各分量 RMSE/MAE/95% 分位数/最大差异。这里的 RMSE 是两份估计之间的差异；各算法相对真值的定位误差由前四个实验计算。结果比较本身不推断两次实验配置相同，使用者应选择对应的算法、参数与数据。

## Colab

使用本目录的 Colab.ipynb，里面有五个实验单元。安装单元把包安装到当前 notebook 的 Python 环境，然后直接 import 同一套示例 run 函数，不再自动建立第二套解释器或通过子进程调用实验。

将 Colab.ipynb 上传到 Colab，选择 github 获取已发布的仓库和分支；它们必须实际包含本次新结构。尚未推送时可以选择 upload，上传 GitHub Download ZIP 或手工压缩的源码 ZIP，压缩包中保留本 README 列出的源码、配置、数据、tests 和安装文件，排除环境、缓存及历史结果。支持 ZIP 根目录即工作目录、包含 kfv_fgo 子目录以及 GitHub 下载的外层目录。也可选择 existing 使用本地目录。重复加载不会删除已有目录。

第五项默认读取 tests/fixtures 下的仿真 MATLAB 参考文件；填写第三项实验生成的 KFV_EKF.npz 路径后即可比较，也可指定自己上传的结果。无需 MATLAB、本地运行时或 MATLAB 适配脚本。最后一个单元将本轮结果 ZIP 保存到 results 下供下载。依赖验证基于 Python 3.12；若当前内核无法安装固定版本的 pyRTKLIB，先选择兼容环境或关闭 INSTALL_REAL。修改或更换代码来源后建议重启内核，避免保留旧模块。

## 验证与发布

```powershell
python -m pytest
```

发布时直接维护并提交本目录的源码、配置、数据、示例、测试、Notebook、安装配置和许可说明，不再维护导出脚本或第二份 Python 副本。GitHub 的源码 ZIP 即可用于分发；dist 不是运行依赖。

重排前源码保存在 .validation/before_simplification.zip，退役文件保存在 .validation/retired-layout，仅供恢复，不参与运行。动态加载、旧通用 CLI 和自动启动 MATLAB 的流程已从维护目录移出。此前验证资料也保留在归档中。

.colab、.validation、results、build、dist、缓存和 egg-info 均忽略版本控制，不随源码发布。本地 Python 环境和验证归档可以继续留在本机。已保存的 MATLAB 参考结果及其来源记录统一保留在 tests/fixtures，运行后新生成的 Python 结果和图放在 results。回归测试覆盖四种 KFV、原模板分支、SWFGO 窗口 1/5/10 及 GNSS/INS 中间量。
