# 项目介绍
本项目通过WiFi信号的信道状态信息（Channel State Information, CSI）进行步态特征提取，进而实现用户识别，大部分算法是对WiFiU[1]的复现，周期提取算法有改动[2]，频谱生成部分参考了[WiHF的开源代码](https://github.com/shulinyang/ReleaseCode)[3]

# 项目使用说明
1. 文件夹前的编号为源码阅读顺序，源码的使用方式在`example.m`中给出
2. `GaitUserID.mlapp`可以通过MATLAB App Designer打开并运行，运行后将出现一个可视化界面，引导用户使用算法实现的各种功能

**请不要随意更改`GaitUserID.mlapp`的位置**，因为在程序启动时，会自动将文件夹下的脚本添加进MATLAB的环境中，如果更改位置可能导致脚本添加环境失败，导致程序无法正常工作

# 源码说明
如果想要阅读源码，可以按照以下顺序进行阅读
- data_loading
- signal_processing
- feature_extraction
- model_training

在`experiments`文件夹中，存放了各种实验功能，以及绘图脚本，模型训练的部分结果存放于`thesis_figures`中


# 运行环境
本项目运行在MATLAB R2022a版本下，并未在R2022a以下的版本测试，因此可能出现sdk版本不兼容的情况，如果出现这种情况，请尝试更新MATLAB版本

# 数据集说明
## 数据要求
本项目在测试时使用[Widar 3.0](http://tns.thss.tsinghua.edu.cn/widar3.0/)中的步态数据进行测试，其文件名格式为`id-a-b-Rx.dat`，其中各字段含义如下：
- 'id': 用户id
- 'a': 路径编号
- 'b': 重复次数，奇数次为朝向发射设备方向，偶数次为远离发射设备方向
- 'Rx': Wi-Fi接收器编号

需要注意的是，由于多普勒频移特征的限制，在收集文件时有特定的要求，**可以处理的`路径-接收器`对为`1-R3`和`3-R6`，且重复次数应为偶数次**

如果是自己收集的数据，只要符合上述要求，也可以直接使用

## 训练数据命名
在使用`GaitUserID.mlapp`提供的模型训练功能时，需要注意将各个用户的数据放在以`userN-XXX`命名的文件夹内，其中`N`为用户编号，其格式为阿拉伯数字，`-XXX`留作备注用，可以是任何符合文件名规范（不包含`.`、`?`等符号）的字符串

# 参考引用
[1] W Wang, A X Liu, M Shahzad. Gait Recognition Using Wifi Signals[C]//Proceedings of the 2016 ACM International Joint Conference on Pervasive and Ubiquitous Computing. Heidelberg Germany: ACM, 2016: 363-373.

[2] https://zhuanlan.zhihu.com/p/21285190

[3] C Li, M Liu, Z Cao. WiHF: Enable User Identified Gesture Recognition with WiFi[C]//IEEE INFOCOM 2020 - IEEE Conference on Computer Communications. Toronto, ON, Canada: IEEE, 2020: 586-595.
