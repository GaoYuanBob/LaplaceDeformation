# LaplaceDeformation
#### 这是基于 MeshLab 2016 的源码，在VS2015上编写的对一个 MeshModel 类型的模型进行 Laplace Deformation 的插件，核心代码只有一个 edit_LaplaceDeformation.h 和一个 edit_LaplaceDeformation.cpp。
使用时只需要在MeshLab界面右侧的模型列表中，用鼠标左键单击要形变的模型（即选为当前模型），即可计算形变，命令行窗口中有输出提示。

需要准备：

1、要形变模型（MeshModel& m）,传入 StartEdit 函数  
2、固定锚点文件 R"(C:\\Users\\Administrator\\Desktop\\fixed_anchor.txt)"，存放固定不动的点的索引，用空格或者换行分隔  
3、移动锚点文件 R"(C:\\Users\\Administrator\\Desktop\\move_anchor_coord.txt)"，每行为 idx x y z 形式，idx 为移动锚点索引，x y z为这个点在形变后的坐标

注：  
1、其中点的索引与 MeshLab 2016 保持一致，从0开始  
2、固定锚点会变为蓝色，移动锚点会变为红色  
3、我有两篇 CSDN 博客详细介绍形变方法和整个流程  
https://blog.csdn.net/bob__yuan/article/details/81778875  
https://blog.csdn.net/bob__yuan/article/details/82191363

下两图为用一个类似于烟囱形状的长方体模型，进行旋转的结果，因为是最基础的 Laplace Deformation 所以`不是保体积的`，旋转多次后会内缩
<div align=center><img src="https://github.com/GaoYuanBob/LaplaceDeformation/raw/master/Sample%20Ressult.png"><div>
<div align=center><img src="https://github.com/GaoYuanBob/LaplaceDeformation/blob/master/Laplace%20Deformation%20-%20Rotation.gif"><div>
