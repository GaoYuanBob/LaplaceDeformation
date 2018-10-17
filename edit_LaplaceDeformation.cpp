#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <meshlab/glarea.h>
#include "edit_LaplaceDeformation.h"
#include <wrap/gl/pick.h>
#include <wrap/qt/gl_label.h>
#include <eigenlib/Eigen/Cholesky>
#include <fstream>

using namespace std;	// MatrixXd 重载了 std::cout << ，直接 cout << m << endl;即可
using namespace vcg;
using Eigen::MatrixXf;
using Eigen::VectorXf;

EditLaplaceDeformationPlugin::EditLaplaceDeformationPlugin()
{
	qFont.setFamily("Helvetica");
	qFont.setPixelSize(12);

	haveToPick = false;
	pickmode = 0; // 0 face 1 vertex
	curFacePtr = 0;
	curVertPtr = 0;
	pIndex = 0;
}

const QString EditLaplaceDeformationPlugin::Info()
{
	return tr("Laplace Deformation on a model!");
}

void EditLaplaceDeformationPlugin::ChimneyRotate(MeshModel &m)
{
	fixed_anchor_idx.clear();
	move_anchor_idx.clear();
	Vertices.clear();
	Faces.clear();
	AdjacentVertices.clear();
	tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
	tri::UpdateSelection<CMeshO>::FaceClear(m.cm);

	int tmp;
	ifstream move_anchor;
	move_anchor.open("C:\\Users\\Administrator\\Desktop\\top_anchor.txt");
	while (move_anchor >> tmp)
		move_anchor_idx.push_back(tmp);
	move_anchor.close();

	ifstream fix_anchor;
	fix_anchor.open("C:\\Users\\Administrator\\Desktop\\bot_anchor.txt");
	while (fix_anchor >> tmp)
		fixed_anchor_idx.push_back(tmp);
	fix_anchor.close();

	// move_anchor设为红色，move_anchor设为蓝色
	//for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	//	if (!vi->IsD())
	//	{
	//		int tmp = vi - m.cm.vert.begin();
	//		if (find(move_anchor_idx.begin(), move_anchor_idx.end(), tmp) != move_anchor_idx.end())
	//			vi->C() = vcg::Color4b(vcg::Color4b::Red);
	//		else if (find(fixed_anchor_idx.begin(), fixed_anchor_idx.end(), tmp) != fixed_anchor_idx.end())
	//			vi->C() = vcg::Color4b(vcg::Color4b::Blue);
	//	}
	
//// test 测试旋转矩阵写得对不对
	//float angle = (45 * M_PI) / 180;	// 每次旋转10度
	//cout << "sin = " << sin(angle) << endl;
	//cout << "cos = " << cos(angle) << endl << endl;
	//for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	//{
	//	int idx = vi - m.cm.vert.begin();
	//	if (find(move_anchor_idx.begin(), move_anchor_idx.end(), idx) != move_anchor_idx.end())
	//	{
	//		float tx = m.cm.vert[idx].P().X(), &ty = m.cm.vert[idx].P().Y();
	//		printf("Before Rotate, (%f, %f)\n", tx, ty);
	//		m.cm.vert[idx].P().X() = tx * cos(angle) + ty * sin(angle);
	//		m.cm.vert[idx].P().Y() = -tx * sin(angle) + ty * cos(angle);
	//		printf("After Rotate, (%f, %f)\n\n", m.cm.vert[idx].P().X(), m.cm.vert[idx].P().Y());
	//	}
	//}
//// test 测试旋转矩阵写得对不对

//// Laplace Deformation - Rotate
	toCaculateAdjacentVertices(&m.cm);
	printf("\nVertices = %d\n", Vertices.size());
	printf("Faces = %d\n\n", Faces.size());

	get_LsTLs_Matrix();
	printf("LsTLs 矩阵计算完成\n");

	int fix_anchors = fixed_anchor_idx.size();
	int move_anchors = move_anchor_idx.size();
	int points_num = Vertices.size();
	VectorXf vx(points_num), vy(points_num);

	for (int i = 0; i < points_num; i++)
	{
		vx[i] = Vertices[i].X();
		vy[i] = Vertices[i].Y();
	}

	bx = L * vx;	// 根据laplace矩阵计算出所有点的的laplace坐标
	by = L * vy;	// 根据laplace矩阵计算出所有点的的laplace坐标

	bx.conservativeResize(points_num + fix_anchors + move_anchors);
	by.conservativeResize(points_num + fix_anchors + move_anchors);


	for (int i = 0; i < fix_anchors; i++)
	{
		bx[i + points_num] = Vertices[fixed_anchor_idx[i]].X();
		by[i + points_num] = Vertices[fixed_anchor_idx[i]].Y();
	}

	float angle = (15.0 * M_PI) / 180.0;	// 每次旋转10度，不能写 90 / 180 * pi，90/180 = 0.。

	// x_new = x * cos + y * sin
	// y_new = -x * sin + y * cos
	for (int i = 0; i < move_anchors; i++)
	{
		float tx = Vertices[move_anchor_idx[i]].X(), ty = Vertices[move_anchor_idx[i]].Y();
		bx[i + points_num + fix_anchors] = tx * cos(angle) + ty * sin(angle);
		by[i + points_num + fix_anchors] = -tx * sin(angle) + ty * cos(angle);
	}
	LsTbx = LsT * bx;
	LsTby = LsT * by;
	printf("LsTb 矩阵计算完成, 通过cholesky分解，解线性方程组 LsTLs * x = LsTb\n\n\n");

	vx_new = LsTLs.llt().solve(LsTbx);
	vy_new = LsTLs.llt().solve(LsTby);

	for (int i = 0; i < m.cm.vert.size(); i++)
	{
		if (find(move_anchor_idx.begin(), move_anchor_idx.end(), i) != move_anchor_idx.end())
		{
			float tx = Vertices[i].X(), ty = Vertices[i].Y();
			m.cm.vert[i].P().X() = tx * cos(angle) + ty * sin(angle);
			m.cm.vert[i].P().Y() = -tx * sin(angle) + ty * cos(angle);
		}
		else
		{
			if (find(fixed_anchor_idx.begin(), fixed_anchor_idx.end(), i) == fixed_anchor_idx.end())	// 既不是move点，又不是固定点
			{
				m.cm.vert[i].P().X() = vx_new[i];
				m.cm.vert[i].P().Y() = vy_new[i];
			}
		}
	}
	//// Laplace Deformation - Rotate

	m.updateDataMask(MeshModel::MM_POLYGONAL);
	tri::UpdateBounding<CMeshO>::Box(m.cm);
	tri::UpdateNormal<CMeshO>::PerVertexNormalizedPerFaceNormalized(m.cm);

	////// 所有update操作
	gla->mvc()->sharedDataContext()->meshInserted(m.id());
	MLRenderingData dt;
	gla->mvc()->sharedDataContext()->getRenderInfoPerMeshView(m.id(), gla->context(), dt);
	suggestedRenderingData(m, dt);
	MLPoliciesStandAloneFunctions::disableRedundatRenderingDataAccordingToPriorities(dt);
	gla->mvc()->sharedDataContext()->setRenderingDataPerMeshView(m.id(), gla->context(), dt);
	gla->update();
	////// 所有update操作
}

// 处理当前选中模型
bool EditLaplaceDeformationPlugin::StartEdit(MeshModel & m, GLArea * _gla, MLSceneGLSharedDataContext* ctx)
{
	gla = _gla;

//// 做烟囱旋转实验 2018.8.31成功，但是不保体积，而且旋转过多面片会有自相交
	//ChimneyRotate(m);
	//return true;
//// 做烟囱旋转实验

	// 为什么模型上来就有选中的点？？必须先ClearS
	/*for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	if (!vi->IsD() && vi->IsS())
	{
	(*vi).C() = vcg::Color4b(vcg::Color4b::Red);
	cout << vi - m.cm.vert.begin() << endl;
	}*/

//// Iniatial 部分
	fixed_anchor_idx.clear();
	move_anchor_idx.clear();
	tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
	tri::UpdateSelection<CMeshO>::FaceClear(m.cm);
	//tri::UpdateFlags<CMeshO>::FaceClearS(m.cm);		// 不对
	//tri::UpdateFlags<CMeshO>::VertexClearS(m.cm);		// 不对

	move_anchor_coord.resize(m.cm.vert.size());
	for (int i = 0; i < move_anchor_coord.size(); i++)
	{
		move_anchor_coord[i].X() = -1;
		move_anchor_coord[i].Y() = -1;
		move_anchor_coord[i].Z() = -1;
	}
//// Iniatial 部分

	m.updateDataMask(MeshModel::MM_ALL);

	// 读入 fixed_anchor.txt 文件，作为固定锚点
	int tmp;
	ifstream fix_a;
	fix_a.open("C:\\Users\\Administrator\\Desktop\\fixed_anchor.txt");
	while (fix_a >> tmp)
		fixed_anchor_idx.push_back(tmp);
	fix_a.close();

	//// 读入 move_anchor.txt 文件，作为移动锚点
	//ifstream move_a;
	//move_a.open("C:\\Users\\Administrator\\Desktop\\move_anchor.txt");
	//while (move_a >> tmp)
	//	move_anchor_idx.push_back(tmp);
	//move_a.close();

	// 读入 move_anchor_coord.txt 文件，作为移动锚点的新坐标 idx x y z

	ifstream move_coord;
	move_coord.open("C:\\Users\\Administrator\\Desktop\\move_anchor_coord.txt");
	while (move_coord >> tmp)
	{
		move_anchor_idx.push_back(tmp);
		move_coord >> move_anchor_coord[tmp].X();
		move_coord >> move_anchor_coord[tmp].Y();
		move_coord >> move_anchor_coord[tmp].Z();
		//if (fabs(m.cm.vert[tmp].P().X() - move_anchor_coord[tmp].X()) > 0.001
		//	&& fabs(m.cm.vert[tmp].P().Y() - move_anchor_coord[tmp].Y()) > 0.001
		//	&& fabs(m.cm.vert[tmp].P().Z() - move_anchor_coord[tmp].Z()) > 0.001)
		//{
		//	printf("m.cm.vert[%d].P() = (%f, %f, %f)\n", tmp, m.cm.vert[tmp].P().X(), m.cm.vert[tmp].P().Y(), m.cm.vert[tmp].P().Z());
		//	printf("move_anchor_coord[%d] = (%f, %f, %f)\n", tmp, move_anchor_coord[tmp].X(), move_anchor_coord[tmp].Y(), move_anchor_coord[tmp].Z());
		//	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
		//}
	}
	move_coord.close();
	
	//int cnt_d = 0;
	//for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	//{
	//	if (vi->IsD())
	//	{
	//		++cnt_d;
	//		printf("点已经被删除\n");
	//	}
	//	else
	//	{
	//		int t = vi - m.cm.vert.begin();
	//		if (find(fixed_anchor_idx.begin(), fixed_anchor_idx.end(), t) != fixed_anchor_idx.end())
	//			vi->C() = vcg::Color4b(vcg::Color4b::Blue);
	//		else if (find(move_anchor_idx.begin(), move_anchor_idx.end(), t) != move_anchor_idx.end())
	//		{
	//			//cout << t << endl;
	//			//cout << m.cm.vert[t].P().Y() << endl;
	//			//cout << move_anchor_coord[t].Y() << endl << endl;

	//			vi->C() = vcg::Color4b(vcg::Color4b::Red);
	//		}
	//	}
	//}
	//printf("\n\n有 %d 个点被删除!\n", cnt_d);
	printf("有 %d 个点被选为固定锚点（蓝色显示）!\n", fixed_anchor_idx.size());
	printf("有 %d 个点被选为移动锚点（红色显示）!\n", move_anchor_idx.size());

	// 进行Laplace形变
	LaplaceDeformation(m);

	// 检查有没有 nan
	//for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); vi++)
	//	printf("%d - (%f, %f, %f)\n", vi - m.cm.vert.begin(), vi->P().X(), vi->P().Y(), vi->P().Z());

////// 所有update操作
	gla->mvc()->sharedDataContext()->meshInserted(m.id());
	MLRenderingData dt;
	gla->mvc()->sharedDataContext()->getRenderInfoPerMeshView(m.id(), gla->context(), dt);
	suggestedRenderingData(m, dt);
	MLPoliciesStandAloneFunctions::disableRedundatRenderingDataAccordingToPriorities(dt);
	gla->mvc()->sharedDataContext()->setRenderingDataPerMeshView(m.id(), gla->context(), dt);
	gla->update();
////// 所有update操作

///////////// 传进来的指针不好使？？？？？// 2018.8.17，上边的代码可以更新模型，但是下边的这段不行
	//ctx->meshInserted(m.id());
	//MLRenderingData dt;
	//ctx->getRenderInfoPerMeshView(m.id(), gla->context(), dt);
	//// 需要自己实现的纯虚函数
	//suggestedRenderingData(m, dt);
	//MLPoliciesStandAloneFunctions::disableRedundatRenderingDataAccordingToPriorities(dt);
	//ctx->setRenderingDataPerMeshView(m.id(), gla->context(), dt);
	//gla->update();
///////////// 传进来的指针不好使？？？？？

	return true;
}

void EditLaplaceDeformationPlugin::EndEdit(MeshModel &/*m*/, GLArea * /*parent*/, MLSceneGLSharedDataContext* /*cont*/)
{
	haveToPick = false;
	pickmode = 0; // 0 face 1 vertex
	curFacePtr = 0;
	curVertPtr = 0;
	pIndex = 0;
}

void EditLaplaceDeformationPlugin::LaplaceDeformation(MeshModel& m)
{
	long t1 = GetTickCount();

	toCaculateAdjacentVertices(&m.cm);
	printf("\nVertices = %d\n", Vertices.size());
	printf("Faces = %d\n\n", Faces.size());

	get_LsTLs_Matrix();
	printf("LsTLs 矩阵计算完成\n");
	get_LsTb_Matrix();
	printf("LsTb 矩阵计算完成, 通过cholesky分解，解线性方程组 LsTLs * x = LsTb\n");

	// 通过cholesky分解，解线性方程组 Ax = b，即 LsTLs * x = LsTb
	// Cholesky分解是把一个对称正定(symmetric, positive definite)的矩阵表示成一个下三角矩阵 L 和其转置 LT 的乘积的分解
	// 需要头文件 #include <Eigen/Cholesky>
	vx_new = LsTLs.ldlt().solve(LsTbx);	// 不知道llt().solve()为什么会出现nan
	//vy_new = LsTLs.ldlt().solve(LsTby);
	//vz_new = LsTLs.ldlt().solve(LsTbz);

	printf("新的坐标计算完成\n");

	setNewCoord(m);
	printf("模型坐标更新完成\n");
	tri::UpdateBounding<CMeshO>::Box(m.cm);
	tri::UpdateNormal<CMeshO>::PerVertexNormalizedPerFaceNormalized(m.cm);

	long t2 = GetTickCount();
	long time_used = (t2 - t1) * 1.0 / 1000;
	int minutes = time_used / 60;
	int seconds = time_used % 60;
	printf("\nLaplace形变完成，用时: %dm %ds\n\n", minutes, seconds);
}

void EditLaplaceDeformationPlugin::toCaculateAdjacentVertices(CMeshO* cm)
{
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;

	for (vi = cm->vert.begin(); vi != cm->vert.end(); ++vi)
		Vertices.push_back(vi->P());

	for (fi = cm->face.begin(); fi != cm->face.end(); ++fi)
		if (!(*fi).IsD())
		{
			std::vector<int> tri;
			for (int k = 0; k < (*fi).VN(); k++)
			{
				int vInd = (int)(fi->V(k) - &*(cm->vert.begin()));
				//if (cm->vert[vInd].IsD())
				//	printf("点已经被删除");
				tri.push_back(vInd);
			}
			Faces.push_back(tri);
		}

	CalculateAdjacentVertices(cm);

	// 输出所有点的近邻数量
	/*for (int i = 0; i < AdjacentVertices.size(); i++)
		printf("第 %d 个点有 %d 个近邻\n", i, AdjacentVertices[i].size());*/

	// 输出所有点周围点的索引
	//for (int i = 0; i < AdjacentVertices.size(); i++)
	//	for (int j = 0; j < AdjacentVertices[i].size(); j++)
	//		printf("AdjacentVertices[%d][%d] = %d\n", i, j, AdjacentVertices[i][j]);
}

void EditLaplaceDeformationPlugin::CalculateAdjacentVertices(CMeshO* cm)
{
	AdjacentVertices.resize(Vertices.size());
	for (size_t i = 0; i < Faces.size(); i++)
	{
		std::vector<int>& t = Faces[i];
		//printf("t0 = %d, t[1] = %, t[2] = %d\n", t[0], t[1], t[2]);
		std::vector<int>& p0list = AdjacentVertices[t[0]];
		std::vector<int>& p1list = AdjacentVertices[t[1]];
		std::vector<int>& p2list = AdjacentVertices[t[2]];

		CMeshO::VertexIterator vi0 = cm->vert.begin() + t[0];
		CMeshO::VertexIterator vi1 = cm->vert.begin() + t[1];
		CMeshO::VertexIterator vi2 = cm->vert.begin() + t[2];
		if (!vi1->IsD() && std::find(p0list.begin(), p0list.end(), t[1]) == p0list.end())
			p0list.push_back(t[1]);
		if (!vi2->IsD() && std::find(p0list.begin(), p0list.end(), t[2]) == p0list.end())
			p0list.push_back(t[2]);
		if (!vi0->IsD() && std::find(p1list.begin(), p1list.end(), t[0]) == p1list.end())
			p1list.push_back(t[0]);
		if (!vi2->IsD() && std::find(p1list.begin(), p1list.end(), t[2]) == p1list.end())
			p1list.push_back(t[2]);
		if (!vi0->IsD() && std::find(p2list.begin(), p2list.end(), t[0]) == p2list.end())
			p2list.push_back(t[0]);
		if (!vi1->IsD() && std::find(p2list.begin(), p2list.end(), t[1]) == p2list.end())
			p2list.push_back(t[1]);
	}
}

void EditLaplaceDeformationPlugin::get_LsTLs_Matrix()
{
	int fix_anchors = fixed_anchor_idx.size();
	int move_anchors = move_anchor_idx.size();
	int points_num = Vertices.size();
	L = MatrixXf(points_num, points_num);
	Ls = MatrixXf(points_num + fix_anchors + move_anchors, points_num);

	for (int i = 0; i < points_num; i++)
	{
		for (int j = 0; j < points_num; j++)
		{
			if (i == j)
			{
				//printf("第 %d 个点有 %d 个近邻\n", i, AdjacentVertices[i].size());
				L(i, j) = AdjacentVertices[i].size();
				Ls(i, j) = AdjacentVertices[i].size();
				continue;
			}

			if (std::find(AdjacentVertices[i].begin(), AdjacentVertices[i].end(), j) != AdjacentVertices[i].end())
			{
				L(i, j) = -1;
				Ls(i, j) = -1;
			}
			else
			{
				L(i, j) = 0;
				Ls(i, j) = 0;
			}
		}
	}

	// 固定锚点
	for (int i = 0; i < fix_anchors; i++)
		for (int j = 0; j < points_num; j++)
		{
			if (j == fixed_anchor_idx[i])
				Ls(points_num + i, j) = 1;
			else
				Ls(points_num + i, j) = 0;
		}

	// 移动锚点
	for (int i = 0; i < move_anchor_idx.size(); i++)
		for (int j = 0; j < points_num; j++)
		{
			if (j == move_anchor_idx[i])
				Ls(points_num + fix_anchors + i, j) = 1;
			else
				Ls(points_num + fix_anchors + i, j) = 0;
		}

	//printf("Laplace矩阵 - L\n");
	//cout << L << endl << endl;

	//printf("添加了锚点信息的Laplace矩阵 - Ls\n");
	//cout << Ls << endl << endl;

	printf("Laplace 矩阵计算完成\n");
	// 矩阵太大了打印出来看不出来什么东西
	// 检查一下
	//for (int i = 0; i < points_num; i++)
	//{
	//	int cnt = 0;
	//	int tmp = 0;
	//	for (int j = 0; j < points_num; j++)
	//	{
	//		if (Ls(i, j) == -1)
	//			++cnt;
	//		else if (Ls(i, j) != 0)
	//			tmp = Ls(i, j);
	//	}
	//	if (cnt != tmp)
	//		printf("Error!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	//}


	LsT = Ls.transpose();
	LsTLs = LsT * Ls;

	// 查看矩阵信息
	//printf("\nLsTLs\n");
	//cout << LsTLs.rows() << "  " << LsTLs.cols() << endl;
	//cout << LsTLs << endl;
	//printf("\n");
}

void EditLaplaceDeformationPlugin::get_LsTb_Matrix()
{
	int fix_anchors = fixed_anchor_idx.size();
	int move_anchors = move_anchor_idx.size();
	int points_num = Vertices.size();
	VectorXf vx(points_num), vy(points_num), vz(points_num);

	for (int i = 0; i < points_num; i++)
	{
		vx[i] = Vertices[i].X();
		vy[i] = Vertices[i].Y();
		vz[i] = Vertices[i].Z();
	}

	// 根据laplace矩阵计算出所有点的的laplace坐标
	bx = L * vx;	// 乘几次是几阶
	by = L * vy;
	bz = L * vz;

	bx.conservativeResize(points_num + fix_anchors + move_anchors);
	by.conservativeResize(points_num + fix_anchors + move_anchors);
	bz.conservativeResize(points_num + fix_anchors + move_anchors);

	// 用形变前坐标对固定锚点坐标进行赋值
	for (int i = 0; i < fix_anchors; i++)
	{
		bx[i + points_num] = Vertices[fixed_anchor_idx[i]].X();
		by[i + points_num] = Vertices[fixed_anchor_idx[i]].Y();
		bz[i + points_num] = Vertices[fixed_anchor_idx[i]].Z();
	}

	// 用形变后坐标对移动锚点坐标进行赋值
	for (int i = 0; i < move_anchors; i++)
	{
		bx[i + points_num + fix_anchors] = move_anchor_coord[move_anchor_idx[i]].X();
		//cout << Vertices[move_anchor_idx[i]].Y() << endl;
		//cout << move_anchor_coord[move_anchor_idx[i]].Y() << endl << endl;
		by[i + points_num + fix_anchors] = move_anchor_coord[move_anchor_idx[i]].Y();
		//cout << Vertices[move_anchor_idx[i]].Y() << endl;
		//cout << move_anchor_coord[move_anchor_idx[i]].Y() << endl << endl;
		bz[i + points_num + fix_anchors] = move_anchor_coord[move_anchor_idx[i]].Z();
	}

	//cout << bx << endl << endl;
	//cout << by << endl << endl;

	// 计算三个轴上的 LsTb 向量
	LsTbx = LsT * bx;
	LsTby = LsT * by;
	LsTbz = LsT * bz;

	//cout << LsTbx << endl << endl;
}

void EditLaplaceDeformationPlugin::setNewCoord(MeshModel& m)
{
	for (int i = 0; i < m.cm.vert.size(); i++)
	{
		if (find(move_anchor_idx.begin(), move_anchor_idx.end(), i) != move_anchor_idx.end())
		{
			m.cm.vert[i].P().X() = move_anchor_coord[i].X();
			//m.cm.vert[i].P().Y() = move_anchor_coord[i].Y();
			//m.cm.vert[i].P().Z() = move_anchor_coord[i].Z();
		}
		else
		{
			if (find(fixed_anchor_idx.begin(), fixed_anchor_idx.end(), i) == fixed_anchor_idx.end())
			{
				m.cm.vert[i].P().X() = vx_new[i];
				//m.cm.vert[i].P().Y() = vy_new[i];
				//m.cm.vert[i].P().Z() = vz_new[i];
			}
		}
	}
}

void EditLaplaceDeformationPlugin::suggestedRenderingData(MeshModel &m, MLRenderingData& dt)
{
	if (m.cm.VN() == 0)
		return;
	MLRenderingData::PRIMITIVE_MODALITY pr = MLRenderingData::PR_SOLID;
	if (m.cm.FN() > 0)
		pr = MLRenderingData::PR_SOLID;

	MLRenderingData::RendAtts atts;
	atts[MLRenderingData::ATT_NAMES::ATT_VERTPOSITION] = true;
	atts[MLRenderingData::ATT_NAMES::ATT_VERTCOLOR] = true;
	atts[MLRenderingData::ATT_NAMES::ATT_VERTNORMAL] = true;
	atts[MLRenderingData::ATT_NAMES::ATT_FACENORMAL] = true;
	atts[MLRenderingData::ATT_NAMES::ATT_FACECOLOR] = true;

	MLPerViewGLOptions opts;
	dt.get(opts);
	opts._sel_enabled = true;
	opts._face_sel = true;
	opts._vertex_sel = true;
	dt.set(opts);

	dt.set(pr, atts);
}

