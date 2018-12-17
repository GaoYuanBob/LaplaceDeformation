#include "edit_LaplaceDeformation_factory.h"
#include "edit_LaplaceDeformation.h"

EditLaplaceDeformationFactory::EditLaplaceDeformationFactory()
{
	editSample = new QAction(QIcon(":/images/icon1_info.png"),"Laplace Deformation", this);
	
	actionList << editSample;
	
	foreach(QAction *editAction, actionList)
		editAction->setCheckable(true); 	
}
	
//gets a list of actions available from this plugin
QList<QAction *> EditLaplaceDeformationFactory::actions() const
{
	return actionList;
}

//get the edit tool for the given action
MeshEditInterface* EditLaplaceDeformationFactory::getMeshEditInterface(QAction *action)
{
	if(action == editSample)
		return new EditLaplaceDeformationPlugin();
	assert(0); //should never be asked for an action that isnt here
	return nullptr;
}

QString EditLaplaceDeformationFactory::getEditToolDescription(QAction *)
{
	return EditLaplaceDeformationPlugin::Info();
}

MESHLAB_PLUGIN_NAME_EXPORTER(EditLaplaceDeformationFactory)
