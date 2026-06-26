# mini-sensor-measurement — 传感器测量与仪器

**从零开始、零依赖的 C 语言实现**，覆盖传感器测量全信号链。

涵盖传感器从物理量到数字量的完整路径：传感原理 → 惠斯通电桥 → 信号调理 →
ADC/DAC 数据转换 → 仪器放大器 → MEMS 传感器 → 测量误差分析（GUM） →
生物传感器电化学 → 智能传感器系统。

基于 Fraden《现代传感器手册》、Pallas-Areny & Webster《传感器与信号调理》、
Doebelin《测量系统》及 IEEE/JCGM 测量标准编写。

## 子模块

| 子模块 | 主题 | 参考课程 |
|--------|------|----------|
| [mini-adc-dac-converter](mini-adc-dac-converter/) | Nyquist 采样、SAR ADC、ΔΣ 调制器、量化理论、动态指标（SNR/SFDR/ENOB） | MIT 6.301, Stanford EE315, Berkeley EE140 |
| [mini-biosensor-chemical](mini-biosensor-chemical/) | 电化学换能机制、Nernst 方程、安培/电位法、葡萄糖生物传感器 | MIT 6.555, UC Berkeley BioE 121 |
| [mini-instrumentation-amp](mini-instrumentation-amp/) | 仪表放大器拓扑（三运放/双运放）、增益结构、CMRR、桥式传感器接口 | Stanford EE247, Berkeley EE105, Michigan EECS 411 |
| [mini-measurement-error](mini-measurement-error/) | 误差分类（偏置/随机/迟滞）、GUM 不确定度传播、校准回归、统计检验 | MIT 6.630, ETH 227-0455, Berkeley EE117 |
| [mini-mems-sensor](mini-mems-sensor/) | MEMS 质量-弹簧-阻尼动力学、Allan 方差、IMU 姿态融合（Mahony/Madgwick/Kalman） | MIT 6.630, Stanford EE247, MIT 16.485 |
| [mini-sensor-principle](mini-sensor-principle/) | 传感器定义、惠斯通电桥分析、传递函数、Johnson/Nyquist 噪声模型 | MIT 6.630, Stanford EE247, Berkeley EE117 |
| [mini-signal-conditioning](mini-signal-conditioning/) | 放大级设计、滤波器、线性化、冷端补偿、4–20 mA 电流环、隔离技术 | MIT 6.003, Stanford EE247, Berkeley EE105 |
| [mini-smart-sensor](mini-smart-sensor/) | 能量采集、Allan 方差拟合、边缘 ML 特征提取（RMS、波峰因子、过零率） | Stanford EE267, ETH 227-0436, MIT 6.630 |

## 设计哲学

1. **从信号到数字** — 每个模块跟踪完整的物理→数字路径：传感元件 → 信号调理 → ADC → 数字处理。
2. **标准驱动** — 严格遵循 IEEE Std 1241（ADC 测试）、IEEE Std 1139/952（噪声/Allan 方差）、JCGM GUM（测量不确定度）、ISO 5725（精度）。
3. **零依赖** — 纯 C99/C11 + `libc` + `libm`，不依赖任何外部库或厂商 SDK。
4. **九层知识体系** — L1 定义 → L2 核心概念 → L3 数学结构 → L4 基本定律 → L5 算法 → L6 经典问题 → L7 工程应用 → L8 进阶专题 → L9 研究前沿。

## 构建方法

```sh
# 构建单个子模块
cd mini-adc-dac-converter && make

# 构建所有子模块
for d in mini-*/; do (cd "$d" && make); done

# 运行测试
for d in mini-*/; do (cd "$d" && make test); done
```

## 项目结构

```
mini-sensor-measurement/
├── mini-adc-dac-converter/        # Nyquist 采样、SAR ADC、ΔΣ 调制器、量化理论、ADC 动态指标
├── mini-biosensor-chemical/       # 电化学换能、Nernst 方程、Butler-Volmer 动力学、葡萄糖传感器
├── mini-instrumentation-amp/      # 三运放/双运放拓扑、增益级分析、桥式传感器接口
├── mini-measurement-error/        # 误差分类学、A类/B类不确定度、校准回归、统计滤波
├── mini-mems-sensor/              # MEMS 动力学、噪声过程分析、六位置标定、IMU 传感器融合
├── mini-sensor-principle/         # 传感器定义、惠斯通电桥、传递函数、噪声模型
├── mini-signal-conditioning/      # 放大级联、滤波器设计、线性化、冷端补偿、4–20 mA 环路
├── mini-smart-sensor/             # 能量采集、Allan 方差拟合、边缘 ML 特征提取
├── .gitignore                     # 构建产物与 IDE 排除项
├── README.md                      # 英文版说明
└── README-CN.md                   # 中文版说明（本文件）
```

## 许可证

MIT — 详见各子模块目录。
