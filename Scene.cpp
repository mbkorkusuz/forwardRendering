#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"


using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
bool visible(double den, double num, double &tE, double &tL);
bool visible(double den, double num, double &tE, double &tL)
{
    
    if (den > 0){ // potentially entering
		double t = num / den;
        if (t > tL){
            return false;
        }
        if (t > tE){
            tE = t;
        }
    }
    else if (den < 0){ // potentially leaving
		double t = num / den;

        if (t < tE){
            return false;
        }
        if (t < tL){
            tL = t;
        }
    }
    else if (num > 0) // line parallel to edge
    {
        return false;
    }
        
    return true;
}

double f01(double x, double y, double &x0, double &x1, double &y0, double &y1)
{
	return x*(y0-y1) + y*(x1-x0) + x0*y1 - x1*y0; 
}
double f12(double x, double y, double &x1, double &x2, double &y1, double &y2)
{
	return x*(y1-y2) + y*(x2-x1) + x1*y2 - x2*y1; 
}
double f20(double x, double y, double &x2, double &x0, double &y2, double &y0)
{
	return x*(y2-y0) + y*(x0-x2) + x2*y0 - x0*y2; 
}

Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

Color Scene::makeBetweenZeroAnd255Color(Color value)
{
	value.r = this->makeBetweenZeroAnd255(value.r);
	value.g = this->makeBetweenZeroAnd255(value.g);
	value.b = this->makeBetweenZeroAnd255(value.b);
	return value;
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "magick " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
	string command2;

	command2 = "rm -f " + ppmFileName;
	system(command2.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	
	vector<Vec3 *> bizimVertices = this->vertices;

	vector<vector<Vec4>> *(verticesForEachMesh) = new vector<vector<Vec4>>();
	(verticesForEachMesh)->resize(this->meshes.size(), vector<Vec4>(bizimVertices.size()));

	vector<vector<double>> *ZDepthVec = new vector<vector<double>>();
	ZDepthVec->resize(camera->verRes, std::vector<double>(camera->horRes));

	for (int i = 0; i < this->meshes.size(); i++){
		for (int j = 0; j < bizimVertices.size(); j++){
			(*verticesForEachMesh)[i][j].x = bizimVertices[j]->x;
			(*verticesForEachMesh)[i][j].y = bizimVertices[j]->y;
			(*verticesForEachMesh)[i][j].z = bizimVertices[j]->z;
			(*verticesForEachMesh)[i][j].t= 1;
			(*verticesForEachMesh)[i][j].colorId = bizimVertices[j]->colorId;
		}	
	}
	
	for (int i = 0; i < camera->verRes; i++){
		for (int j = 0; j < camera->horRes; j++){
			(*ZDepthVec)[i][j] = INT_MAX;
		}
	}
	

	Matrix4 Mcam;
	Matrix4 Mloc;
	Mloc.values[0][0] = 1; Mloc.values[1][1] = 1; Mloc.values[2][2] = 1; Mloc.values[3][3] = 1;
	Mloc.values[0][3] = -camera->position.x ; Mloc.values[1][3] = -camera->position.y; Mloc.values[2][3] = -camera->position.z;
	Matrix4 Muvw;
	Muvw.values[3][3] = 1;
	Muvw.values[0][0] = camera->u.x; Muvw.values[0][1] = camera->u.y; Muvw.values[0][2] = camera->u.z; 
	Muvw.values[1][0] = camera->v.x; Muvw.values[1][1] = camera->v.y; Muvw.values[1][2] = camera->v.z; 
	Muvw.values[2][0] = camera->w.x; Muvw.values[2][1] = camera->w.y; Muvw.values[2][2] = camera->w.z; 
	Mcam = multiplyMatrixWithMatrix(Muvw, Mloc);

	Matrix4 Mper;
	if (camera->projectionType == 0){
		// orthographic
		Mper.values[0][0] = 2/(camera->right - camera->left);
		Mper.values[0][3] = -(camera->right + camera->left)/(camera->right - camera->left);
		Mper.values[1][1] = 2/(camera->top - camera->bottom);
		Mper.values[1][3] = -(camera->top + camera->bottom)/(camera->top - camera->bottom);
		Mper.values[2][2] = -2/(camera->far - camera->near);
		Mper.values[2][3] =  -(camera->far + camera->near)/(camera->far - camera->near);
		Mper.values[3][3] = 1;
	}
	else {
		// pers
		Mper.values[0][0] = 2*(camera->near)/(camera->right - camera->left);
		Mper.values[0][2] = (camera->right + camera->left)/(camera->right - camera->left);
		Mper.values[1][1] = 2*(camera->near)/(camera->top - camera->bottom);
		Mper.values[1][2] = (camera->top + camera->bottom)/(camera->top - camera->bottom);
		Mper.values[2][2] = -(camera->far + camera->near)/(camera->far - camera->near);
		Mper.values[2][3] = (-2*(camera->far * camera->near))/(camera->far - camera->near);
		Mper.values[3][2] = -1;
	}

	Matrix4 MVP;
	MVP.values[0][0] = float((camera->horRes)) / 2; 
	MVP.values[0][3] = float((camera->horRes - 1))  / 2; 
	MVP.values[1][1] = float((camera->verRes)) / 2; 
	MVP.values[1][3] = float((camera->verRes - 1))  / 2; 
	MVP.values[2][2] = 0.5; 
	MVP.values[2][3] = 0.5; 
	
	
	for (int i = 0; i < this->meshes.size(); i++){
		Mesh currentMesh = *(this->meshes[i]);
		Vec3 u, v, w;
		double te, tl, dx, dy, dz, num;
		bool vsb;

		// identity matrix olusturuldu
		Matrix4 modelingMatrix;
		modelingMatrix.values[0][0] = 1;
		modelingMatrix.values[1][1] = 1;
		modelingMatrix.values[2][2] = 1;
		modelingMatrix.values[3][3] = 1;
		
		Matrix4 transformationMatrix;
		transformationMatrix.values[0][0] = 1;
		transformationMatrix.values[1][1] = 1;
		transformationMatrix.values[2][2] = 1;
		transformationMatrix.values[3][3] = 1;

		Matrix4 rxMatrix;
		rxMatrix.values[0][0] = 1;
		rxMatrix.values[1][1] = 1;
		rxMatrix.values[2][2] = 1;
		rxMatrix.values[3][3] = 1;

		Matrix4 M;
		M.values[3][3] = 1;

		Matrix4 MEksibir;
		MEksibir.values[3][3] = 1;


		for (int j = 0; j < currentMesh.numberOfTransformations; j++){
			
			char transformationTypeOfMesh = currentMesh.transformationTypes[j];
			int transformationIdOfMesh = currentMesh.transformationIds[j];

			if (transformationTypeOfMesh == 't'){
				// translation matrix
				transformationMatrix.values[0][3] = this->translations[transformationIdOfMesh-1]->tx;
				transformationMatrix.values[1][3] = this->translations[transformationIdOfMesh-1]->ty;
				transformationMatrix.values[2][3] = this->translations[transformationIdOfMesh-1]->tz;
			}
			else if (transformationTypeOfMesh == 's'){
				// scaling matrix
				transformationMatrix.values[0][0] = this->scalings[transformationIdOfMesh-1]->sx;
				transformationMatrix.values[1][1] = this->scalings[transformationIdOfMesh-1]->sy;
				transformationMatrix.values[2][2] = this->scalings[transformationIdOfMesh-1]->sz;
			}
			else if (transformationTypeOfMesh == 'r'){
				// rotation matrix
				u.x = this->rotations[transformationIdOfMesh-1]->ux;
				u.y = this->rotations[transformationIdOfMesh-1]->uy;
				u.z = this->rotations[transformationIdOfMesh-1]->uz;

				if (abs(u.x) <= abs(u.y) && abs(u.x) <= abs(u.z)){
					v.x = 0;
					v.y = -u.z;
					v.z = u.y;
				}
				else if(abs(u.y) <= abs(u.x) && abs(u.y) <= abs(u.z)){
					v.y = 0;
					v.x = -u.z;
					v.z = u.x;
				}
				else{
					v.z = 0;
					v.x = -u.y;
					v.y = u.x;
				}
				v = normalizeVec3(v);
				w = crossProductVec3(u,v);
				w = normalizeVec3(w);

				rxMatrix.values[1][1] = cos((this->rotations[transformationIdOfMesh-1]->angle * M_PI) / 180.0);
				rxMatrix.values[1][2] = -sin((this->rotations[transformationIdOfMesh-1]->angle * M_PI) / 180.0);
				rxMatrix.values[2][1] = sin((this->rotations[transformationIdOfMesh-1]->angle * M_PI) / 180.0);
				rxMatrix.values[2][2] = cos((this->rotations[transformationIdOfMesh-1]->angle * M_PI) / 180.0);

				M.values[0][0] = u.x; M.values[0][1] = u.y; M.values[0][2] = u.z;
				M.values[1][0] = v.x; M.values[1][1] = v.y; M.values[1][2] = v.z;
				M.values[2][0] = w.x; M.values[2][1] = w.y; M.values[2][2] = w.z;

				MEksibir.values[0][0] = u.x; MEksibir.values[1][0] = u.y; MEksibir.values[2][0] = u.z;
				MEksibir.values[0][1] = v.x; MEksibir.values[1][1] = v.y; MEksibir.values[2][1] = v.z;
				MEksibir.values[0][2] = w.x; MEksibir.values[1][2] = w.y; MEksibir.values[2][2] = w.z;

				

				transformationMatrix = multiplyMatrixWithMatrix(MEksibir,multiplyMatrixWithMatrix(rxMatrix, M));
			}
			
			modelingMatrix = multiplyMatrixWithMatrix(transformationMatrix, modelingMatrix);
			transformationMatrix = Matrix4(getIdentityMatrix());
			
		}
		
		for (int j = 0; j < bizimVertices.size(); j++){
			
			(*verticesForEachMesh)[i][j] = multiplyMatrixWithVec4(modelingMatrix, (*verticesForEachMesh)[i][j]);
			
			(*verticesForEachMesh)[i][j] = multiplyMatrixWithVec4(Mcam, (*verticesForEachMesh)[i][j]);
			
			(*verticesForEachMesh)[i][j] = multiplyMatrixWithVec4(Mper, (*verticesForEachMesh)[i][j]);
			
		}
 
		for (int j = 0; j < bizimVertices.size(); j++){
			(*verticesForEachMesh)[i][j].x /= (*verticesForEachMesh)[i][j].t;
			(*verticesForEachMesh)[i][j].y /= (*verticesForEachMesh)[i][j].t;
			(*verticesForEachMesh)[i][j].z /= (*verticesForEachMesh)[i][j].t;
			(*verticesForEachMesh)[i][j].t /= (*verticesForEachMesh)[i][j].t;
		}
		// perspective division tamamlandi

		vector<bool> trianglesbool;

		for(int j = 0;  j < currentMesh.triangles.size(); j++){
			Vec3 a,b,c;
			Vec3 normal;
			Vec3 backfaceV;
			Vec3 CameraPos;

			CameraPos.x = 0; CameraPos.y = 0; CameraPos.z = 0;

			a.x = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].x; a.y = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].y; a.z = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].z;
			b.x = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].x; b.y = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].y; b.z = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].z;
			c.x = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].x; c.y = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].y; c.z = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].z;

			backfaceV = normalizeVec3(subtractVec3(CameraPos, a));
			normal = crossProductVec3(subtractVec3(c, b), subtractVec3(a,b));

			if (dotProductVec3(normal, backfaceV) > 0) {
				trianglesbool.push_back(false);
			}
			else{
				trianglesbool.push_back(true);
			}
		}

		

		std::vector<Vec4> linesX0;
		std::vector<Vec4> linesX1;

		for (int j = 0; j < currentMesh.numberOfTriangles; j++)
		{	
			Vec4 X0((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1]);
			Vec4 X1((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1]);
			
			linesX0.push_back(X0);
			linesX1.push_back(X1);
			
			X0 = Vec4((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1]);
			X1 = Vec4((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1]);

			linesX0.push_back(X0);
			linesX1.push_back(X1);

			X0 = Vec4((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1]);
			X1 = Vec4((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1]);

			linesX0.push_back(X0);
			linesX1.push_back(X1);
			
		}
		
		
		// clipping
		if(currentMesh.type == 0){
		
		for (int j = 0; j < linesX0.size(); j++){
			// xmin = -1, xmax = 1 ... diye varsaydik  
			te = 0;
			tl = 1;
			vsb = false;
			dx = linesX1[j].x - linesX0[j].x;
			dy = linesX1[j].y - linesX0[j].y;
			dz = linesX1[j].z - linesX0[j].z;

			if (visible(dx, -1 - linesX0[j].x, te, tl)){
				if (visible(-dx, -1 + linesX0[j].x, te, tl)){
					if (visible(dy, -1 - linesX0[j].y, te, tl)){
						if (visible(-dy, -1 + linesX0[j].y, te, tl)){
							if (visible(dz, -1 - linesX0[j].z, te, tl)){
								if (visible(-dz, -1 + linesX0[j].z, te, tl)){
									vsb = true;

									if (tl < 1){
										linesX1[j].x = linesX0[j].x + dx*tl;
										linesX1[j].y = linesX0[j].y + dy*tl;
										linesX1[j].z = linesX0[j].z + dz*tl;
									}
									if (te > 0){
										linesX0[j].x = linesX0[j].x + dx*te;
										linesX0[j].y = linesX0[j].y + dy*te;
										linesX0[j].z = linesX0[j].z + dz*te;
									}
									
								}
							}
						}
					}
				}
			}
		}
	}	
		
		for (int j = 0; j < bizimVertices.size(); j++){
			(*verticesForEachMesh)[i][j] = multiplyMatrixWithVec4(MVP, (*verticesForEachMesh)[i][j]);
		}

		for (int j = 0; j < linesX0.size(); j++)
		{
			linesX0[j] = multiplyMatrixWithVec4(MVP, linesX0[j]);
			linesX1[j] = multiplyMatrixWithVec4(MVP, linesX1[j]);
		}
		
	
		// viewport transformation tamamlandi

		//rasterization
		
		double x0, x1, y0, y1, z0, z1, x0new, y0new, x1new, y1new;
		double d, slop, zDep;
		int y, x, p;
		p = -1;
		
		Color c, dc;
		if(currentMesh.type == 0){
			// line rasterization
			for (int j = 0; j < linesX0.size(); j++){

				if (j%3 == 0){
					p++;
				}
				if (trianglesbool[p] == false && this->cullingEnabled == true){
					continue;
				}
				if(linesX0[j].x < -0.5 || linesX0[j].y < -0.5 || linesX1[j].x < -0.5 || linesX1[j].y < -0.5){
					continue;
				}
				if(linesX0[j].x > camera->horRes || linesX0[j].y > camera->verRes || linesX1[j].x > camera->horRes || linesX1[j].y > camera->verRes){
					continue;
				}
			
				// a-b 
				if (linesX0[j].x < linesX1[j].x){
					
					x0 = linesX0[j].x;
					x1 = linesX1[j].x;
					y0 = linesX0[j].y;
					y1 = linesX1[j].y;
					z0 = linesX0[j].z;
					z1 = linesX1[j].z;
					if(j%3 == 0){
						x0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].x;
						x1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].x;
						y0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].y;
						y1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].y;
					}
					if(j%3 == 1){
						x0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].x;
						x1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].x;
						y0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].y;
						y1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].y;
					}
					if(j%3 == 2){
						x0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].x;
						x1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].x;
						y0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].y;
						y1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].y;
					}

					if (x1-x0 == 0){
						
						if (y1-y0 > 0){
							slop = INT_MAX;
						}
						else if (y1 - y0 < 0){
							slop = INT_MIN;
						}
						else{
							slop = 0;
						}
					}
					else{
						
						slop = ((double)y1-y0)/((double)x1-x0);
					}
					c = *(this->colorsOfVertices[linesX0[j].colorId-1]);
					if(x0new < x0){
						double scrics;
						scrics = (x0 - x0new)/(x1new - x0new);
						c.r += scrics * ((this->colorsOfVertices[linesX1[j].colorId-1])->r - c.r);
						c.g += scrics * ((this->colorsOfVertices[linesX1[j].colorId-1])->g - c.g);
						c.b += scrics * ((this->colorsOfVertices[linesX1[j].colorId-1])->b - c.b);

					}

					if(slop < -1 || slop > 1){
						if(y1-y0 != 0){
							
						dc.r = ((this->colorsOfVertices[linesX1[j].colorId-1])->r - (this->colorsOfVertices[linesX0[j].colorId-1])->r) / (y1new - y0new);
						dc.g = ((this->colorsOfVertices[linesX1[j].colorId-1])->g - (this->colorsOfVertices[linesX0[j].colorId-1])->g) / (y1new - y0new);
						dc.b = ((this->colorsOfVertices[linesX1[j].colorId-1])->b - (this->colorsOfVertices[linesX0[j].colorId-1])->b) / (y1new - y0new);
						
						}
					}
					else{
						if(x1-x0 != 0){
						dc.r = ((this->colorsOfVertices[linesX1[j].colorId-1])->r - (this->colorsOfVertices[linesX0[j].colorId-1])->r) / (x1new - x0new);
						dc.g = ((this->colorsOfVertices[linesX1[j].colorId-1])->g - (this->colorsOfVertices[linesX0[j].colorId-1])->g) / (x1new - x0new);
						dc.b = ((this->colorsOfVertices[linesX1[j].colorId-1])->b - (this->colorsOfVertices[linesX0[j].colorId-1])->b) / (x1new - x0new);
						}
					}
				}
				else{
					x1 = linesX0[j].x;
					x0 = linesX1[j].x;
					y1 = linesX0[j].y;
					y0 = linesX1[j].y;
					z1 = linesX0[j].z;
					z0 = linesX1[j].z;
					

					if(j%3 == 0){
						x1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].x;
						x0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].x;
						y1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].y;
						y0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].y;
					}
					if(j%3 == 1){
						x1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].x;
						x0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].x;
						y1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[1] - 1].y;
						y0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].y;
					}
					if(j%3 == 2){
						x1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].x;
						x0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].x;
						y1new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[2] - 1].y;
						y0new = (*verticesForEachMesh)[i][currentMesh.triangles[p].vertexIds[0] - 1].y;
					}
					if (x1-x0 == 0){
						if (y1-y0 > 0){
							slop = INT_MAX;
						}
						else if (y1 - y0 < 0){
							slop = INT_MIN;
						}
						else{
							slop = 0;
						}
					}
					else{
						slop = ((double)y1-y0)/((double)x1-x0);
					}
					c = *(this->colorsOfVertices[linesX1[j].colorId-1]);
					if(x0new < x0){
						double scrics;
						scrics = (x0 - x0new)/(x1new - x0new);
						c.r += scrics * ((this->colorsOfVertices[linesX0[j].colorId-1])->r - c.r);
						c.g += scrics * ((this->colorsOfVertices[linesX0[j].colorId-1])->g - c.g);
						c.b += scrics * ((this->colorsOfVertices[linesX0[j].colorId-1])->b - c.b);

					}
					if(slop < -1 || slop > 1){
						if(y1-y0 != 0){
						dc.r = ((this->colorsOfVertices[linesX0[j].colorId-1])->r - (this->colorsOfVertices[linesX1[j].colorId-1])->r) / (y1new - y0new);
						dc.g = ((this->colorsOfVertices[linesX0[j].colorId-1])->g - (this->colorsOfVertices[linesX1[j].colorId-1])->g) / (y1new - y0new);
						dc.b = ((this->colorsOfVertices[linesX0[j].colorId-1])->b - (this->colorsOfVertices[linesX1[j].colorId-1])->b) / (y1new - y0new);
						}
					}
					else{
						if(x1-x0 != 0){
						dc.r = ((this->colorsOfVertices[linesX0[j].colorId-1])->r - (this->colorsOfVertices[linesX1[j].colorId-1])->r) / (x1new - x0new);
						dc.g = ((this->colorsOfVertices[linesX0[j].colorId-1])->g - (this->colorsOfVertices[linesX1[j].colorId-1])->g) / (x1new - x0new);
						dc.b = ((this->colorsOfVertices[linesX0[j].colorId-1])->b - (this->colorsOfVertices[linesX1[j].colorId-1])->b) / (x1new - x0new);
						}
					}
				}
				zDep = z0;
					
				if(slop > 1){
					//for y0 to y1
					d = (x0 - x1) + (0.5 * (y1 - y0));
					x = x0;
					for (int k = y0; k <= y1; k++)
					{
						
						
						if ((*ZDepthVec)[x][k] > zDep){
							(*ZDepthVec)[x][k] = zDep;
							
							this->assignColorToPixel(x,k,makeBetweenZeroAnd255Color(c));
						}
						
						zDep += (z1-z0)/(y1-y0);
						if (d < 0){
							x += 1;
							d += (x0 - x1) + (y1 - y0);
						}
						else{
							d += x0 - x1;
						}
						c.r += dc.r;
						c.g += dc.g;
						c.b += dc.b;
					}
					
				}
				else if(slop>=0){
					
					//for x0 to x1
					d = (y0 - y1) + (0.5 * (x1 - x0));
					y = y0;
					for (int k = x0; k <= x1; k++)
					{
						if ((*ZDepthVec)[k][y] > zDep){
							(*ZDepthVec)[k][y] = zDep;
							this->assignColorToPixel(k, y, makeBetweenZeroAnd255Color(c));
						}
						zDep += (z1-z0)/(x1-x0);
						if (d < 0){
							y += 1;
							d += (y0 - y1) + (x1 - x0);
						}
						else{
							d += y0 - y1;
						}
						c.r += dc.r;
						c.g += dc.g;
						c.b += dc.b;
					}
				}
				else if(slop>=-1){
					//for x0 to x1
					// y yi azalt
					d = (y1 - y0) + (0.5 * (x1 - x0));
					y = y0;

					for (int k = x0; k <= x1; k++)
					{
						if ((*ZDepthVec)[k][y] > zDep){
							(*ZDepthVec)[k][y] = zDep;
							this->assignColorToPixel(k, y, makeBetweenZeroAnd255Color(c));
						}
						zDep += (z1-z0)/(x1-x0);
						
						if (d < 0){
							y -= 1;
							d += (y1 - y0) + (x1 - x0);
						}
						else{
							d += y1 - y0;
						}
						c.r += dc.r;
						c.g += dc.g;
						c.b += dc.b;
					}
				}
				else{
					// for y0 to y1 eksilerek gidecek
					// xleri arttir
					d = (x0 - x1) + (0.5 * (y0 - y1));
					x = x0;
					
					for (int k = y0; k >= y1; k--)
					{
						if ((*ZDepthVec)[x][k] > zDep){
							(*ZDepthVec)[x][k] = zDep;
							this->assignColorToPixel(x,k,makeBetweenZeroAnd255Color(c));
						}
						zDep += (z1-z0)/(y0-y1);
						
						if (d < 0){
							x += 1;
							d += (x0 - x1) + (y0 - y1);
						}
						else{
							d += x0 - x1;
						}
						c.r -= dc.r;
						c.g -= dc.g;
						c.b -= dc.b;
					}
				}
			}
		}
		else{
			
			// triangle rasterization
			int xmin, ymin, xmax, ymax;
			double alpha, beta, theta, zdep;
			Color cc;

			double X0, X1, X2, Y0, Y1, Y2;
			
			for (int j = 0; j < currentMesh.numberOfTriangles; j++){

				
				if (trianglesbool[j] == false && this->cullingEnabled == true){
					continue;
				}
				
				X0 = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].x;
				X1 = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].x;
				X2 = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].x;
				Y0 = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].y;
				Y1 = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].y;
				Y2 = (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].y;

				xmin = min((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].x, min((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].x, (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].x));
				xmax = max((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].x, max((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].x, (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].x));
				ymin = min((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].y, min((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].y, (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].y));
				ymax = max((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].y, max((*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].y, (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].y));

				cc = *(this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].colorId-1]);
				for (int ty = ymin; ty <= ymax; ty++){
					for (int tx = xmin; tx <= xmax; tx++){
						if (tx < 0 || ty < 0 || tx > camera->horRes || ty > camera->verRes){
							continue;
						}
						if(tx == camera->horRes || ty == camera->verRes){
							continue;
						}

						alpha = f12(tx, ty, X1, X2, Y1, Y2) / f12(X0, Y0, X1, X2, Y1, Y2);
						beta = f20(tx, ty, X2, X0, Y2, Y0) / f20(X1, Y1, X2, X0, Y2, Y0);
						theta = f01(tx, ty, X0, X1, Y0, Y1) / f01(X2, Y2, X0, X1, Y0, Y1);

						zdep = alpha * (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].z + beta * (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].z + theta * (*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].z;
						
						if ((*ZDepthVec)[ty][tx] < zdep){
							continue;
						}
						
						if (alpha >= 0 && beta >= 0 && theta >= 0){
							(*ZDepthVec)[ty][tx] = zdep;
							cc.r = alpha * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].colorId - 1])->r + beta * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].colorId - 1])->r + theta * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].colorId - 1])->r;
							cc.g = alpha * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].colorId - 1])->g + beta * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].colorId - 1])->g + theta * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].colorId - 1])->g;
							cc.b = alpha * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[0] - 1].colorId - 1])->b + beta * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[1] - 1].colorId - 1])->b + theta * (this->colorsOfVertices[(*verticesForEachMesh)[i][currentMesh.triangles[j].vertexIds[2] - 1].colorId - 1])->b;
							this->assignColorToPixel(tx, ty, makeBetweenZeroAnd255Color(cc));
						}	
					}	
				}
			}	
		}
	}
}


