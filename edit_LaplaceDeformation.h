#pragma once

#include <vector>
#include <QObject>
#include <common/interfaces.h>
#include <vcg/space/deprecated_point3.h>

#include <eigenlib/Eigen/Dense>

class EditLaplaceDeformationPlugin : public QObject, public MeshEditInterface
{
	Q_OBJECT
	Q_INTERFACES(MeshEditInterface)
		
public:
	Eigen::MatrixXf L, Ls, LsT, LsTLs, LsTbx, LsTby, LsTbz;	// L: 原始laplace矩阵（方阵），LS: 加入了锚点（point_num + anchors, point_num）大小， T：转置
	Eigen::VectorXf bx, by, bz;
	Eigen::VectorXf vx_new, vy_new, vz_new;

	std::vector<int> fixed_anchor_idx, move_anchor_idx;	// 锚点分为固定点和移动点
	std::vector<vcg::Point3f> move_anchor_coord;		// 形变后移动锚点位置
	std::vector<vcg::Point3f> Vertices;
	std::vector<std::vector<int>> Faces;
	std::vector<std::vector<int>> AdjacentVertices;

	EditLaplaceDeformationPlugin();
    virtual ~EditLaplaceDeformationPlugin() = default;	// resharper C++

	static const QString Info();

    bool StartEdit(MeshModel &/*m*/, GLArea * /*parent*/, MLSceneGLSharedDataContext* /*cont*/) override;	// resharper C++ 加了个 "override"
    void EndEdit(MeshModel &/*m*/, GLArea * /*parent*/, MLSceneGLSharedDataContext* /*cont*/) override;
	void Decorate(MeshModel &/*m*/, GLArea * /*parent*/, QPainter *p) override {};
	void mousePressEvent(QMouseEvent *, MeshModel &, GLArea * ) override {};
	void mouseMoveEvent(QMouseEvent *, MeshModel &, GLArea * ) override {};
	void mouseReleaseEvent(QMouseEvent *event, MeshModel &/*m*/, GLArea *) override {};
	void keyReleaseEvent(QKeyEvent *, MeshModel &, GLArea *) override {};

	void LaplaceDeformation(MeshModel&);	// 依次调用下边函数

	void toCaculateAdjacentVertices(CMeshO* cm);
	void CalculateAdjacentVertices(CMeshO* cm);

	void get_LsTLs_Matrix();	// 获取laplace矩阵以及LTLILT矩阵
	void get_LsTb_Matrix();		// 获取 b 矩阵，x' = LTLILT * b

	void setNewCoord(MeshModel&);	// 更新模型坐标，anchor点单独按照形变计算

	void suggestedRenderingData(MeshModel &, MLRenderingData& dt) override;
	
	void ChimneyRotate(MeshModel &);

private:
    QPoint cur;
	QFont qFont;
    bool haveToPick;
	int pickmode;
	CMeshO::FacePointer   curFacePtr;
	CMeshO::VertexPointer curVertPtr;
	int pIndex;
	GLArea * gla;
};
