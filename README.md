NIfTI 磁共振体数据转 DICOM 用户手册
1. 概述
nifti_to_dicom.py 将心脏幅度 NIfTI 体数据（单个 4D 文件，或一组有序的 3D 文件）转换为经典单帧 MR DICOM（MRImageStorage），适用于在不支持 4D/增强 MR 的查看器或后处理软件中加载。

核心特性

按命令行顺序严格读取输入，不做配准、重采样或插值

支持单相位导出（无需时序信息）与多相位导出（需提供 RR 间期或触发时间）

默认将平面内像素重排为“前部朝上”（anterior-up）显示布局；可用 --native-layout 保留原始体素顺序

所有相位共享一套全局 16 位线性量化（RescaleSlope/RescaleIntercept），保证相位间强度可比

输出目录必须为空，绝不覆盖已有文件

生成 conversion_report.json 记录几何、量化、时序与 UID 信息

不适用范围：速度/相位图、DWI、分割掩膜、增强 MR（Enhanced MR）导出。这些场景请使用专门的转换工具。

2. 系统要求
Python 3.8+

依赖包：

bash
pip install numpy nibabel "pydicom>=3.0"
注意：pydicom 必须为 3.0 或更高版本（脚本使用了 enforce_file_format=True 等新接口）。

3. 快速开始
3.1 单相位导出（最简单）
bash
python nifti_to_dicom.py \
  --input cardiac_4d.nii.gz \
  --reference patient_ref.dcm \
  --output ./dicom_out \
  --single-phase
默认导出第 1 个相位；--single-phase 3 导出第 3 个相位

单相位导出不需要时序参数

3.2 多相位导出（已知 RR 间期）
bash
python nifti_to_dicom.py \
  --input cardiac_4d.nii.gz \
  --reference patient_ref.dcm \
  --output ./dicom_out \
  --rr-ms 800 \
  --start-ms 0
假设相位在 RR 间期内均匀采样：

text
times[t] = start_ms + t * rr_ms / N   (t = 0 … N-1)
3.3 多相位导出（已知每个相位的触发时间）
bash
python nifti_to_dicom.py \
  --input phase01.nii.gz phase02.nii.gz phase03.nii.gz \
  --reference patient_ref.dcm \
  --output ./dicom_out \
  --trigger-times-ms 0,150,320,480,650
多个 3D 文件的相位顺序严格按命令行顺序，请务必核对。

4. 命令行参数详解
参数	必填	默认值	说明
--input	是	—	一个 4D NIfTI 文件，或一组有序 3D NIfTI 文件（空格分隔）
--reference	是	—	同患者的参考 DICOM，用于继承患者/研究信息
--output	是	—	输出目录；必须不存在或为空
--trigger-times-ms	否	—	逗号分隔的 R 波后延迟（ms），每个相位一个；与 --rr-ms 互斥
--rr-ms	否	—	标称 RR 间期（ms）；假设相位均匀采样；与 --trigger-times-ms 互斥
--single-phase [N]	否	不启用	仅导出一个相位，1-based；不带值时默认 1；此时时序参数可选
--start-ms	否	0	第一个相位的延迟，仅与 --rr-ms 配合使用
--organization	否	single-series	single-series：所有相位放入同一 Series；per-phase：每个相位一个 Series
--description	否	Cardiac NIfTI derived	SeriesDescription（最长 64 字符）
--series-number	否	900	基准 SeriesNumber；per-phase 时逐相位递增
--same-frame-of-reference	否	关闭	复用参考 DICOM 的 FrameOfReferenceUID；仅当 NIfTI 保留参考患者坐标时使用
--assume-mm	否	关闭	将 NIfTI 中未知空间单位显式视为 mm；使用前务必确认
--copy-mr-parameters	否	关闭	从参考 DICOM 复制采集参数（TR/TE/翻转角等）；仅当参考与源采集匹配时使用
--native-layout	否	关闭	保持原始平面内体素顺序，不做“前部朝上”重排
4.1 时序参数规则
单相位：--rr-ms 和 --trigger-times-ms 都可省略；若提供 --trigger-times-ms，长度可以是 1（即所选相位时间）或输入相位总数（脚本会自动选取对应项）。

多相位：必须提供 --rr-ms 或 --trigger-times-ms 之一，否则报错。

--trigger-times-ms 必须是有限、非负、严格递增的数值。

--rr-ms 必须有限且为正；且 start_ms + (N-1)*rr_ms/N < rr_ms，否则报错。
--start-ms 非零但未提供 --rr-ms 时会报错。

4.2 布局参数说明
默认（anterior-up）：脚本根据 NIfTI 仿射的 LPS 方向，自动选择前后（AP）投影最大的平面内轴作为行轴，并翻转使前部朝上。

若前后方向垂直于切片平面（例如矢状面重建），脚本会报错，提示使用 --native-layout 或先重采样到轴位/斜位平面。

--native-layout：保留原始 (i,j) 体素顺序，不做翻转或转置。

5. 使用场景与示例
5.1 单个 4D 文件，导出全部相位（单 Series）
bash
python nifti_to_dicom.py \
  --input cines.nii.gz \
  --reference ref.dcm \
  --output out_all \
  --rr-ms 850 \
  --start-ms 20 \
  --organization single-series
5.2 多个 3D 文件，每相位独立 Series
bash
python nifti_to_dicom.py \
  --input p1.nii.gz p2.nii.gz p3.nii.gz p4.nii.gz \
  --reference ref.dcm \
  --output out_per_phase \
  --trigger-times-ms 0,200,400,600 \
  --organization per-phase \
  --series-number 1000 \
  --description "Cine derived"
5.3 只导出第 2 个相位，且不提供时序
bash
python nifti_to_dicom.py \
  --input cines.nii.gz \
  --reference ref.dcm \
  --output out_phase2 \
  --single-phase 2
5.4 保留原始布局并复用 FrameOfReferenceUID
bash
python nifti_to_dicom.py \
  --input cines.nii.gz \
  --reference ref.dcm \
  --output out_native \
  --single-phase 1 \
  --native-layout \
  --same-frame-of-reference
5.5 复制参考 MR 采集参数
bash
python nifti_to_dicom.py \
  --input cines.nii.gz \
  --reference ref.dcm \
  --output out_params \
  --rr-ms 900 \
  --copy-mr-parameters
--copy-mr-parameters 仅复制 ScanningSequence、RepetitionTime、EchoTime、FlipAngle、MagneticFieldStrength 等字段；请确认参考 DICOM 与 NIfTI 来源为同一采集，否则会引入错误元数据。

6. 输出文件结构
text
<output>/
├── phase_001/
│   ├── IM_00001.dcm
│   ├── IM_00002.dcm
│   └── ...
├── phase_002/
│   └── ...
└── conversion_report.json
每个相位一个子目录（即使 single-series 组织方式也按相位分目录存放文件；Series 归属由 DICOM 内的 SeriesInstanceUID 决定）

每个切片一个 DICOM 文件，文件名 IM_<slice:05d>.dcm

conversion_report.json 包含：

输入形状、相位数、输出文件数

组织方式、触发时间、单相位选择信息

时序假设（unknown / uniform RR sampling / explicit user times）

输入 LPS 仿射、输出 LPS 仿射、输出形状

输出到输入的索引映射矩阵（output_to_input_indices）

anterior_up 标志

slope / intercept / max_quantization_error

生成的 Series UID 列表

frame_of_reference_reused 标志

7. 处理流程与算法说明
几何读取：Nibabel 仿射（RAS）转换为 DICOM 患者坐标（LPS），即 diag(-1,-1,1,1) @ affine。要求 sform_code 或 qform_code 非零。

单位处理：支持 mm、meter、micron；未知单位必须加 --assume-mm。

仿射校验：拒绝含剪切的仿射（要求方向矩阵正交）；要求所有相位共享相同形状与仿射。

布局重排：根据前后方向选择行轴与列轴，生成索引映射矩阵 mapping；像素按该映射转置/翻转，无插值。

全局量化：

扫描所有相位得到全局 [lo, hi]

slope = (hi - lo) / 65535（若 hi == lo 则 slope = 1）

intercept = lo

存储值 = round((value - intercept) / slope)，裁剪到 [0, 65535]

DICOM 写出：

从参考 DICOM 复制患者/研究级标签（PatientName、StudyDate、AccessionNumber 等）

设置 ImageType = ['DERIVED','SECONDARY','OTHER']

多相位时设置 TemporalPositionIdentifier、NumberOfTemporalPositions，若触发时间等间隔则设置 TemporalResolution

提供时序时设置 TriggerTime；提供 --rr-ms 时设置 NominalInterval

PositionReferenceIndicator 置空

SliceThickness / SpacingBetweenSlices 使用重建体素宽度（非恢复的采集层厚）

8. 限制与不支持的情况
❌ 速度/相位图、DWI、分割掩膜、增强 MR

❌ 含剪切的仿射（需先重采样到正交网格）

❌ 前后方向垂直于切片平面且未加 --native-layout（需重采样到轴位/斜位）

❌ 行/列超过 65535（经典 DICOM 限制）

❌ 含 NaN/Inf 的数据（需先清理）

❌ 多相位输入形状或物理几何不一致（需先配准/重采样）

❌ 输出目录非空（脚本绝不覆盖已有文件）

9. 常见错误与排查
错误信息	原因与处理
Output directory must be empty; existing output is never overwritten.	输出目录已存在且非空；更换目录或手动清空
Use one 4D file or a list of 3D files.	--input 混用了 4D 与 3D，或单个文件不是 4D
For a single 3D volume or singleton 4D input, add --single-phase.	单相位输入需显式加 --single-phase
Multi-phase export requires --rr-ms or --trigger-times-ms.	多相位导出必须提供时序参数
--start-ms is only used with --rr-ms.	--start-ms 不能单独使用
Trigger times must be finite, nonnegative and increasing...	检查 --trigger-times-ms 的数值与数量
Last phase must be within RR interval; check --start-ms.	start_ms + (N-1)*rr_ms/N 必须小于 rr_ms
NIfTI has no coded sform/qform; patient geometry is unknown.	NIfTI 缺少空间信息；需重新导出带 sform/qform 的文件
Unknown spatial units: use --assume-mm only after checking.	确认单位后加 --assume-mm
Sheared affine: resample to an orthogonal voxel grid first.	仿射含剪切；先重采样
All phases must share shape and physical geometry.	相位间形状或仿射不一致；先配准/重采样
NaN/Inf present; clean the data explicitly before conversion.	数据含非有限值；先清理
Anterior direction is perpendicular to the slice plane...	加 --native-layout 或先重采样到轴位/斜位
Reference has no FrameOfReferenceUID.	参考 DICOM 缺少 FrameOfReferenceUID；换参考或去掉 --same-frame-of-reference

10. 转换后验证建议

脚本本身在结束时提示：

Done. Re-import and verify geometry, intensities, and temporal grouping in target viewer.

建议在目标查看器中核对：

几何：切片位置、方向、像素间距是否与源数据一致（对照 conversion_report.json 中的 output_affine_lps_mm）

强度：反量化后 stored * slope + intercept 是否与源 NIfTI 一致（最大误差约 slope/2）

时序分组：single-series 下是否能按 TemporalPositionIdentifier 正确分组；per-phase 下各 Series 是否独立

显示方向：默认前部朝上；查看器可能覆盖存储方向，以实际显示为准

11. 最小工作示例（完整流程）
bash
# 1. 安装依赖
pip install numpy nibabel "pydicom>=3.0"

# 2. 单相位快速验证
python nifti_to_dicom.py \
  --input cines.nii.gz \
  --reference ref.dcm \
  --output ./test_single \
  --single-phase

# 3. 检查输出
ls test_single/phase_001 | head
cat test_single/conversion_report.json

# 4. 多相位正式导出
python nifti_to_dicom.py \
  --input cines.nii.gz \
  --reference ref.dcm \
  --output ./final \
  --rr-ms 800 \
  --start-ms 0 \
  --organization single-series
本手册基于 nifti_to_dicom.py 当前版本编写；若脚本参数或行为有更新，请以脚本内 --help 输出与源码注释为准。