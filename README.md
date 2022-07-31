# turbulenceResponseModel

## 简介

## 安装方法

1. 下载分支
```shell
cd $FOAM_SRC/phaseSystemModels/twoPhaseEuler/phaseCompressibleTurbulenceModels
git clone https://github.com/fightingxiaoxiao/turbulenceResponseModel.git
```

2. 添加头文件
在`phaseCompressibleTurbulenceModels.C`中附加头文件：
```cpp
// -------------------------------------------------------------------------- //
// Additional models
// -------------------------------------------------------------------------- //

#include "kineticTheoryModel.H"
makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, kineticTheoryModel);

#include "phasePressureModel.H"
makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, phasePressureModel);

// 新增部分
#include "turbulenceResponseModel.H" 
makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, turbulenceResponseModel);

// ************************************************************************* //
```

3. 修改构建脚本
`Make/files`修改为：
```shell
phaseCompressibleTurbulenceModels.C
phasePressureModel/phasePressureModel.C
turbulenceResponseModel/turbulenceResponseModel.C # 增加这行

kineticTheoryModels/kineticTheoryModel/kineticTheoryModel.C
...
```
