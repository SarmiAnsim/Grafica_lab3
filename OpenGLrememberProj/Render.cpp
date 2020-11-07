#include "Render.h"

#include <sstream>
#include <iostream>
#include <iomanip>

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include "GUItextRectangle.h"


#include <chrono>
#include <Windows.h>
#include <cmath>
#define _USE_MATH_DEFINES 
#include <math.h>	
#include <vector>
#define Vector std::vector


class Point
{
public: Point(double A, double B, double C)
{
	x = A;
	y = B;
	z = C;
}
	  Point(double* pnt)
	  {
		  x = pnt[0];
		  y = pnt[1];
		  z = pnt[2];
	  }

	  Point(Vector<double> pnt)
	  {
		  x = pnt[0];
		  y = pnt[1];
		  z = pnt[2];
	  }

	  Point operator- (Point p) { return { p.x - x, p.y - y, p.z - z }; }

	  void NewValue(double* pnt)
	  {
		  x = pnt[0];
		  y = pnt[1];
		  z = pnt[2];
	  }

	  double Abs()
	  {
		  return sqrt(x * x + y * y + z * z);
	  }

	  double Abs(int O)
	  {
		  if (O == 0)
			  return sqrt(y * y + z * z);
		  if (O == 1)
			  return sqrt(x * x + z * z);
		  if (O == 2)
			  return sqrt(y * y + x * x);
	  }

	  double scal(Point p)
	  {
		  return p.x * x + p.y * y + p.z * z;
	  }

	  Point vect(Point p)
	  {
		  return { y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x };
	  }

	  Point rottX(double angle)
	  {
		  return { x, y * cos(angle) + z * sin(angle), -y * sin(angle) + z * cos(angle) };
	  }

	  Point rottZ(double angle)
	  {
		  return { x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle), z };
	  }

	  double AX()
	  {
		  if (z >= 0 && y >= 0)
			  return M_PI / 2 + asin(z / Abs(0)); // Назад, нижний сектор
		  if (z >= 0 && y < 0)
			  return M_PI / 2 + asin(z / Abs(0)); // Назад, верхний сектор
		  if (z < 0 && y < 0)
			  return M_PI / 2 - asin(abs(z) / Abs(0)); // Вперед, нижний сектор
		  if (z < 0 && y >= 0)
			  return M_PI / 2 - asin(abs(z) / Abs(0)); // Вперед, верхний сектор
		  return M_PI / 2;
	  }
	  double AZ()
	  {
		  if (y >= 0 && x >= 0)
			  return asin(y / Abs(2)) + M_PI / 2; // Вперед, верхний сектор
		  if (y >= 0 && x < 0)
			  return -M_PI / 2 - asin(y / Abs(2)); // Назад, нижний сектор
		  if (y < 0 && x < 0)
			  return asin(abs(y) / Abs(2)) - M_PI / 2; // Назад, верхний сектор
		  if (y < 0 && x >= 0)
			  return M_PI / 2 - asin(abs(y) / Abs(2)); // Вперед, нижний сектор
		  return 0;
	  }

	  void norm()
	  {
		  double l = sqrt(x * x + y * y + z * z);
		  if (l)
		  {
			  x /= l;
			  y /= l;
			  z /= l;
		  }
	  }

	  double x, y, z;
};

void DrowRocket(Point A, Point angle);
void first_task(double delta_time);
void second_task(double detail);

Point Norm(double* A, double* B, double* C);

bool changePoint = true;

//Point Norm(double A[3], double B[3], double C[3]);
//
//double text_o(int k, int i, int side);

bool textureMode = true;
bool lightMode = true;

bool show_norm = false;

bool rocket_stup = false;

//класс для настройки камеры
class CustomCamera : public Camera
{
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	
	//значния масеры по умолчанию
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera()
	{
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //создаем объект камеры


//Класс для настройки света
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	
	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //создаем источник света

Vector<Vector<Vector<double>>> Pnts = {
	{{0,0,0}, {0,1,1}, {0,2,1}, {0,3,0}},
	{{1,0,1}, {1,1,1}, {1,2,1}, {1,3,1}},
	{{2,0,1}, {2,1,1}, {2,2,1}, {2,3,1}},
	{{3,0,0}, {3,1,1}, {3,2,1}, {3,3,0}}
};

Vector<Vector<Vector<double>>> rez;
Vector<Vector<Vector<Vector<double>>>> norm;

//старые координаты мыши
int mouseX = 0, mouseY = 0;
int I = -1, J = -1;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой левой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	if (!OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y; 

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double diffX = 5, diffY = 5;
		if(I == -1)
			for(int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
				{
					double k = 0, x = 0, y = 0;
					double z = Pnts[i][j][2];

					if (r.direction.Z() == 0)
						k = 0;
					else
						k = (z - r.origin.Z()) / r.direction.Z();

					x = k * r.direction.X() + r.origin.X();
					y = k * r.direction.Y() + r.origin.Y();


					if ((diffX > abs(Pnts[i][j][0] - x) && abs(Pnts[i][j][0] - x) <= 0.2) && (diffY > abs(Pnts[i][j][1] - y) && abs(Pnts[i][j][1] - y) <= 0.2))
					{
						diffX = abs(Pnts[i][j][0] - x);
						diffY = abs(Pnts[i][j][1] - y);
						I = i;
						J = j;
					}
				}
		if(I != -1)
			Pnts[I][J][2] += 0.02 * dy;
		if(!OpenGL::isKeyPressed('S'))
			changePoint = !changePoint;
	}

	if (OpenGL::isKeyPressed('S') && !OpenGL::isKeyPressed(VK_LBUTTON))
		changePoint = !changePoint;

	if (!OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		I = -1;
		J = -1;
	}
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

double detail = 1;

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}

	if (key == 'N')
	{
		show_norm = !show_norm;
	}

	if (key == 'K')
	{
		Pnts.clear();
		Pnts = {
		{{0,0,0}, {0,1,1}, {0,2,1}, {0,3,0}},
		{{1,0,1}, {1,1,1}, {1,2,1}, {1,3,1}},
		{{2,0,1}, {2,1,1}, {2,2,1}, {2,3,1}},
		{{3,0,0}, {3,1,1}, {3,2,1}, {3,3,0}}
		};
		detail = 1;
		changePoint = !changePoint;
	}

	if (key == '0')
	{
		++detail;
		changePoint = !changePoint;
	}

	if (key == '9' && detail >= 2)
	{
		--detail;
		changePoint = !changePoint;
	}

	if (key == 'Z')
	{
		rocket_stup = !rocket_stup;
	}
}

void keyUpEvent(OpenGL *ogl, int key)
{
	
}

GLuint texId[3];

//выполняется перед первым рендером
void initRender(OpenGL* ogl)
{
	//настройка текстур

	//4 байта на хранение пикселя
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//настройка режима наложения текстур
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//включаем текстуры
	glEnable(GL_TEXTURE_2D);


	//массив трехбайтных элементов (R G B)
	RGBTRIPLE* texarray;

	//массив символов, (высота*ширина*4 4, потомучто выше, мы указали использовать по 4 байта на пиксель текстуры - R G B A)
	char* texCharArray;
	int texW, texH;
	OpenGL::LoadBMP("texture1.bmp", &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	RGBTRIPLE* texarray2;
	char* texCharArray2;
	int texW2, texH2;
	OpenGL::LoadBMP("texture2.bmp", &texW2, &texH2, &texarray2);
	OpenGL::RGBtoChar(texarray2, texW2, texH2, &texCharArray2);

	RGBTRIPLE* texarray3;
	char* texCharArray3;
	int texW3, texH3;
	OpenGL::LoadBMP("texture3.bmp", &texW3, &texH3, &texarray3);
	OpenGL::RGBtoChar(texarray3, texW3, texH3, &texCharArray3);

	//генерируем ИД для текстуры
	glGenTextures(3, texId);
	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId[0]);
	//загружаем текстуру в видеопямять, в оперативке нам больше она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);
	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId[1]);
	//загружаем текстуру в видеопямять, в оперативке нам больше она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW2, texH2, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray2);
	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId[2]);
	//загружаем текстуру в видеопямять, в оперативке нам больше она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW3, texH3, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray3);
	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	//отчистка памяти
	free(texCharArray);
	free(texarray);
	free(texCharArray2);
	free(texarray2);
	free(texCharArray3);
	free(texarray3);

	//OpenGL::LoadBMP("texture1.bmp", &texW, &texH, &texarray);
	//OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	////генерируем ИД для текстуры
	//glGenTextures(1, );
	////биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	//glBindTexture(GL_TEXTURE_2D, texID2);

	////загружаем текстуру в видеопямять, в оперативке нам больше она не нужна
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

	////отчистка памяти
	//free(texCharArray);
	//free(texarray);

	////наводим шмон
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH);


	// задать параметры освещения
	// параметр GL_LIGHT_MODEL_TWO_SIDE -
	// 0 - лицевые и изнаночные рисуются одинаково(по умолчанию),
	// 1 - лицевые и изнаночные обрабатываются разными режимами
	// соответственно лицевым и изнаночным свойствам материалов.
	// параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение,
	// не зависящее от сточников
	// по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}

auto end_render = std::chrono::steady_clock::now();

void Render(OpenGL *ogl)
{
	auto cur_time = std::chrono::steady_clock::now();
	auto deltatime = cur_time - end_render;
	double delta = 1.0 * std::chrono::duration_cast<std::chrono::microseconds>(deltatime).count() / 1000000;
	end_render = cur_time;

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);


	//альфаналожение
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//настройка материала
	GLfloat amb[] = { 0.2, 0.2, 0.2, 1. };
	GLfloat dif[] = { 0.5, 0.5, 0.5, 1. };
	GLfloat spec[] = { 0.6, 0.6, 0.6, 1. };
	GLfloat sh = 0.1f * 256;


	//фоновая
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//дифузная
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//зеркальная
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//размер блика
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);
	//===================================
	//Прогать тут  

	first_task(delta);
	second_task(detail);



   //Сообщение вверху экрана

	
	glMatrixMode(GL_PROJECTION);	//Делаем активной матрицу проекций. 
	                                //(всек матричные операции, будут ее видоизменять.)
	glPushMatrix();   //сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
	glLoadIdentity();	  //Загружаем единичную матрицу
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //врубаем режим ортогональной проекции

	glMatrixMode(GL_MODELVIEW);		//переключаемся на модел-вью матрицу
	glPushMatrix();			  //сохраняем текущую матрицу в стек (положение камеры, фактически)
	glLoadIdentity();		  //сбрасываем ее в дефолт

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //классик моего авторства для удобной работы с рендером текста.
	rec.setSize(300, 230);
	rec.setPosition(10, ogl->getHeight() - 230 - 10);


	std::stringstream ss;
	ss << "T - вкл/выкл текстур" << std::endl;
	ss << "L - вкл/выкл освещение" << std::endl;
	ss << "F - Свет из камеры" << std::endl;
	ss << "G - двигать свет по горизонтали" << std::endl;
	ss << "G+ЛКМ двигать свет по вертекали" << std::endl;
	ss << "Детализация поверхности: 9(+)/0(-)" << std::endl;
	ss << "S+ЛКМ - двигать точки без перерисовки" << std::endl;
	ss << "K - восстановить точки и детализацию" << std::endl;
	ss << "N - показать/скрыть нормали (" << norm.size()*norm[0].size()*2 << "шт)" << std::endl;
	ss << "Z - ракету на взлет/запустить" << std::endl;
	ss << "Коорд. света: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "Коорд. камеры: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "Параметры камеры: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //восстанавливаем матрицы проекции и модел-вью обратьно из стека.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}

double f1(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * t * t * t - 3 * t * t + 1) + p4 * (-2 * t * t * t + 3 * t * t) + r1 * (t * t * t - 2 * t * t + t) + r4 * (t * t * t - t * t);
}

double f2(double p1, double p2, double p3, double p4, double t)
{
	return (1 - t) * (1 - t) * (1 - t) * p1 + 3 * t * (1 - t) * (1 - t) * p2 + 3 * t * t * (1 - t) * p3 + t * t * t * p4;
}

double t_max = 0;


long double fact(int N)
{
	if (N < 0)
		return 0;
	if (N == 0)
		return 1;
	else
		return N * fact(N - 1);
}

double polyn_Bernstein(int N, int i, double t)
{
	return fact(N) / (fact(i) * fact(N - i)) * pow(t, i) * pow(1 - t, N - i);
}

double f3(int N, int M, double u, double v, std::vector<std::vector<std::vector<double>>> P, int cord)
{
	double S = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			S += polyn_Bernstein(N - 1, i, u) * polyn_Bernstein(M - 1, j, v) * P[i][j][cord];
		}
	}
	return S;
}


double text(double coord, int ax)
{
	return (coord - rez[0][0][ax]) / (rez[rez.size() - 1][rez[0].size() - 1][ax] - rez[0][0][ax]);
}

void second_task(double detail)
{

	Vector<Vector<double>> rez1;
	Vector<double> rez2;

	Vector<Vector<Vector<double>>> norm1;
	Vector<Vector<double>> norm2;
	Vector<double> norm3;

	if (changePoint)
	{
		norm.clear();
		rez.clear();
		for (double u = 0; u <= 1.0001; u += 0.1/ detail)
		{
			rez1.clear();
			for (double v = 0; v <= 1.0001; v += 0.1/ detail)
			{
				rez2.clear();
				double P[3] = { 0 };

				P[0] = f3(4, 4, u, v, Pnts, 0);
				P[1] = f3(4, 4, u, v, Pnts, 1);
				P[2] = f3(4, 4, u, v, Pnts, 2);

				rez2.push_back(P[0]);
				rez2.push_back(P[1]);
				rez2.push_back(P[2]);
				rez1.push_back(rez2);

			}
			rez.push_back(rez1);
		}

		for (int i = 0; i < rez.size() - 1; i++)
		{
			norm1.clear();
			for (int j = 0; j < rez[i].size() - 1; j++)
			{
				double P1[3] = { rez[i][j][0], rez[i][j][1], rez[i][j][2] };
				double P2[3] = { rez[i + 1][j][0], rez[i + 1][j][1], rez[i + 1][j][2] };
				double P3[3] = { rez[i][j + 1][0], rez[i][j + 1][1], rez[i][j + 1][2] };
				double P4[3] = { rez[i + 1][j + 1][0], rez[i + 1][j + 1][1], rez[i + 1][j + 1][2] };

				norm2.clear();

				norm3.clear();
				Point P = Norm(P1, P2, P3);
				norm3.push_back(P.x);
				norm3.push_back(P.y);
				norm3.push_back(P.z);
				norm2.push_back(norm3);

				norm3.clear();
				P = Norm(P2, P4, P3);
				norm3.push_back(P.x);
				norm3.push_back(P.y);
				norm3.push_back(P.z);
				norm2.push_back(norm3);

				norm1.push_back(norm2);
			}
			norm.push_back(norm1);
		}

		changePoint = !changePoint;
	}

	glBindTexture(GL_TEXTURE_2D, texId[0]);
	glPointSize(5);
	glColor3d(0.5, 0.5, 0.5);
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < rez.size() - 1; i++)
	{ 
		for (int j = 0; j < rez[i].size() - 1; j++)
		{
			double P1[3] = { rez[i][j][0], rez[i][j][1], rez[i][j][2] };
			double P2[3] = { rez[i + 1][j][0], rez[i + 1][j][1], rez[i + 1][j][2] };
			double P3[3] = { rez[i][j + 1][0], rez[i][j + 1][1], rez[i][j + 1][2] };
			double P4[3] = { rez[i + 1][j + 1][0], rez[i + 1][j + 1][1], rez[i + 1][j + 1][2] };

			glNormal3d(norm[i][j][0][0], norm[i][j][0][1], norm[i][j][0][2]);
			glTexCoord2d(text(P1[0], 0), text(P1[1], 1));
			glVertex3dv(P1);
			glTexCoord2d(text(P2[0], 0), text(P2[1], 1));
			glVertex3dv(P2);
			glTexCoord2d(text(P3[0], 0), text(P3[1], 1));
			glVertex3dv(P3);

			glNormal3d(norm[i][j][1][0], norm[i][j][1][1], norm[i][j][1][2]);
			glTexCoord2d(text(P2[0], 0), text(P2[1], 1));
			glVertex3dv(P2);
			glTexCoord2d(text(P3[0], 0), text(P3[1], 1));
			glVertex3dv(P3);
			glTexCoord2d(text(P4[0], 0), text(P4[1], 1));
			glVertex3dv(P4);

		}
	}
	glEnd();

	if (show_norm)
	{
		glDisable(GL_LIGHTING);
		glColor3d(1, 0, 1);
		glBegin(GL_LINES);
		for (int i = 0; i < rez.size() - 1; i++)
		{
			for (int j = 0; j < rez[i].size() - 1; j++)
			{
				double P[3] = { 0 };
				double P1[3] = { rez[i][j][0], rez[i][j][1], rez[i][j][2] };
				double P2[3] = { rez[i + 1][j][0], rez[i + 1][j][1], rez[i + 1][j][2] };
				double P3[3] = { rez[i][j + 1][0], rez[i][j + 1][1], rez[i][j + 1][2] };
				double P4[3] = { rez[i + 1][j + 1][0], rez[i + 1][j + 1][1], rez[i + 1][j + 1][2] };

				glNormal3d(norm[i][j][0][0], norm[i][j][0][1], norm[i][j][0][2]);
				glVertex3dv(P1);
				P[0] = P1[0] + norm[i][j][0][0] / 2;
				P[1] = P1[1] + norm[i][j][0][1] / 2;
				P[2] = P1[2] + norm[i][j][0][2] / 2;
				glVertex3dv(P);
				glVertex3dv(P2);
				P[0] = P2[0] + norm[i][j][0][0] / 2;
				P[1] = P2[1] + norm[i][j][0][1] / 2;
				P[2] = P2[2] + norm[i][j][0][2] / 2;
				glVertex3dv(P);
				glVertex3dv(P3);
				P[0] = P3[0] + norm[i][j][0][0] / 2;
				P[1] = P3[1] + norm[i][j][0][1] / 2;
				P[2] = P3[2] + norm[i][j][0][2] / 2;
				glVertex3dv(P);

				glNormal3d(norm[i][j][1][0], norm[i][j][1][1], norm[i][j][1][2]);
				glVertex3dv(P2);
				P[0] = P2[0] + norm[i][j][1][0] / 2;
				P[1] = P2[1] + norm[i][j][1][1] / 2;
				P[2] = P2[2] + norm[i][j][1][2] / 2;
				glVertex3dv(P);
				glVertex3dv(P3);
				P[0] = P3[0] + norm[i][j][1][0] / 2;
				P[1] = P3[1] + norm[i][j][1][1] / 2;
				P[2] = P3[2] + norm[i][j][1][2] / 2;
				glVertex3dv(P);
				glVertex3dv(P4);
				P[0] = P4[0] + norm[i][j][1][0] / 2;
				P[1] = P4[1] + norm[i][j][1][1] / 2;
				P[2] = P4[2] + norm[i][j][1][2] / 2;
				glVertex3dv(P);

			}
		}
		glEnd();
	}

	glPushAttrib(GL_DEPTH_BUFFER_BIT); // сохраняем стейт глубины
	glDepthFunc(GL_ALWAYS);

	glDisable(GL_LIGHTING);
	/*glPointSize(10);*/
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			double P1[3] = { Pnts[i][j][0], Pnts[i][j][1], Pnts[i][j][2] };

			if (i != 3)
			{
				double P2[3] = { Pnts[i+1][j][0], Pnts[i+1][j][1], Pnts[i+1][j][2] };

				glVertex3dv(P1);
				glVertex3dv(P2);
			}

			if (j != 3)
			{
				double P3[3] = { Pnts[i][j+1][0], Pnts[i][j+1][1], Pnts[i][j+1][2] };

				glVertex3dv(P1);
				glVertex3dv(P3);
			}
		}
	glEnd();


	glPointSize(10);
	glColor3d(1, 0, 0);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			double P[3] = { Pnts[i][j][0], Pnts[i][j][1], Pnts[i][j][2] };
			glVertex3dv(P);
		}
	glEnd();
	glEnable(GL_TEXTURE_2D);

	glPopAttrib(); // восстанавливаем стайт
}

void first_task(double delta_time)
{

	t_max += delta_time / 5;
	if (t_max > 1) t_max = 0;

	double Ermit[2][6][3] = {
		{ { 24,17,30 }, { 26, 31, 35}, { 15, 38, 37}, { 29,36,27 },{0,0,0},{0,0,0} },
		{ { 25,37,31 }, { 26, 34, 22}, { 34, 24, 38}, { 21,34,27 },{0,0,0},{0,0,0} }
	};
	double Beze[2][4][3] = {
		{ { 0,0,0 }, { -5, -5, 25}, { -50, -30, 16}, { -25,-75,42 } },
		{ { 5,7,9 }, { 6, 4, 8}, { 6, 7, 7}, { 4,7,9 } }
		// { 0,0,0 }, { 0.6, -0.8, 20}, { -9, 13, 22}, { -15,-10,25 } 
	};

	for (int i = 0; i < 2; i++)
	{
		Ermit[i][4][0] = Ermit[i][1][0] - Ermit[i][0][0];
		Ermit[i][4][1] = Ermit[i][1][1] - Ermit[i][0][1];
		Ermit[i][4][2] = Ermit[i][1][2] - Ermit[i][0][2];

		Ermit[i][5][0] = Ermit[i][2][0] - Ermit[i][3][0];
		Ermit[i][5][1] = Ermit[i][2][1] - Ermit[i][3][1];
		Ermit[i][5][2] = Ermit[i][2][2] - Ermit[i][3][2];
	}

	double P1[3], P2[3];

	glDisable(GL_LIGHTING);
	glLineWidth(3);
	for (int i = 0; i < 2; i++)
	{
		glColor3d(0, 1, 0);
		glBegin(GL_LINE_STRIP);
		for (double t = 0; t <= 1.0001; t += 0.01)
		{
			P1[0] = f1(Ermit[i][0][0], Ermit[i][3][0], Ermit[i][4][0], Ermit[i][5][0], t);
			P1[1] = f1(Ermit[i][0][1], Ermit[i][3][1], Ermit[i][4][1], Ermit[i][5][1], t);
			P1[2] = f1(Ermit[i][0][2], Ermit[i][3][2], Ermit[i][4][2], Ermit[i][5][2], t);
			glVertex3dv(P1);
		}
		glEnd();

		glColor3d(0, 0, 1);
		glBegin(GL_LINE_STRIP);
		for (double t = 0; t <= 1.0001; t += 0.01)
		{
			P2[0] = f2(Beze[i][0][0], Beze[i][1][0], Beze[i][2][0], Beze[i][3][0], t);
			P2[1] = f2(Beze[i][0][1], Beze[i][1][1], Beze[i][2][1], Beze[i][3][1], t);
			P2[2] = f2(Beze[i][0][2], Beze[i][1][2], Beze[i][2][2], Beze[i][3][2], t);
			glVertex3dv(P2);
		}
		glEnd();
	}


	glLineWidth(1);
	for (int i = 0; i < 2; i++)
	{
		glColor3d(1, 0, 1);
		glBegin(GL_LINES);
		glVertex3dv(Ermit[i][0]);
		glVertex3dv(Ermit[i][1]);

		glVertex3dv(Ermit[i][2]);
		glVertex3dv(Ermit[i][3]);

		glEnd();

		glColor3d(0, 1, 1);
		glBegin(GL_LINE_STRIP);
		glVertex3dv(Beze[i][0]);
		glVertex3dv(Beze[i][1]);
		glVertex3dv(Beze[i][2]);
		glVertex3dv(Beze[i][3]);
		glEnd();

		glPointSize(10);
		glColor3d(1, 0, 0);
		glBegin(GL_POINTS);
		glVertex3dv(Ermit[i][0]);
		glVertex3dv(Ermit[i][1]);
		glVertex3dv(Ermit[i][2]);
		glVertex3dv(Ermit[i][3]);

		glVertex3dv(Beze[i][0]);
		glVertex3dv(Beze[i][1]);
		glVertex3dv(Beze[i][2]);
		glVertex3dv(Beze[i][3]);
		glEnd();
	}

	Point A0(0, 0, 0), A1(0, 0, 0);
	Point D(0, 0, 0), G(0, 0, 0);
	for (double t = 0; t <= 2 * t_max; t += 0.005)
	{
		A0 = A1;
		if (t <= 1)
		{
			P2[0] = f2(Beze[0][0][0], Beze[0][1][0], Beze[0][2][0], Beze[0][3][0], t);
			P2[1] = f2(Beze[0][0][1], Beze[0][1][1], Beze[0][2][1], Beze[0][3][1], t);
			P2[2] = f2(Beze[0][0][2], Beze[0][1][2], Beze[0][2][2], Beze[0][3][2], t);
		}
		else
		{
			--t;
			P2[0] = f2(Beze[0][3][0], Beze[0][2][0], Beze[0][1][0], Beze[0][0][0], t);
			P2[1] = f2(Beze[0][3][1], Beze[0][2][1], Beze[0][1][1], Beze[0][0][1], t);
			P2[2] = f2(Beze[0][3][2], Beze[0][2][2], Beze[0][1][2], Beze[0][0][2], t);
			++t;
		}
		A1.NewValue(P2);
		D = A1 - A0;
		D.norm();
		G = { D.AX(), 0, D.AZ() };
	}

	glEnable(GL_LIGHTING);
	if(!rocket_stup)
		DrowRocket(P2, G);
	else
		DrowRocket({ 0,0,0 }, { 0,0,0 });
}

double changefire = 0;
Vector<Vector<double>> track;

void DrowFire(Vector<Vector<double>> points, Point A, Point angle)
{
	Vector <Vector<Vector<double>>> fire = {
		{{points[0][0] + 0.1, points[0][1]-0.1, points[0][2]}, 
		{points[0][0], points[0][1] - 0.2, points[0][2]}, 
		{points[1][0], points[1][1] + 0.2, points[1][2]},
		{points[1][0] + 0.1, points[1][1] + 0.1, points[1][2]}},

		{{points[0][0] + 0.2, points[0][1], points[0][2]},
		{points[0][0] + 0.2, points[0][1] - 0.2, -3 - changefire},
		{points[1][0] + 0.2, points[1][1] + 0.2, -3 - changefire},
		{points[1][0] + 0.2, points[1][1], points[1][2]}},

		{{points[2][0] - 0.2, points[2][1], points[2][2]},
		{points[2][0] - 0.2, points[2][1] - 0.2, -3 - changefire},
		{points[3][0] - 0.2, points[3][1] + 0.2, -3 - changefire},
		{points[3][0] - 0.2, points[3][1], points[3][2]}},

		{{points[2][0] - 0.1, points[2][1] - 0.1, points[2][2]},
		{points[2][0], points[2][1] - 0.2, points[2][2]},
		{points[3][0], points[3][1] + 0.2, points[3][2]},
		{points[3][0] - 0.1, points[3][1] + 0.1, points[3][2]}}
	};

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Point P = fire[i][j];
			P = P.rottX(angle.x);
			P = P.rottZ(angle.z);
			fire[i][j][0] = P.x + A.x;
			fire[i][j][1] = P.y + A.y;
			fire[i][j][2] = P.z + A.z;
		}
	}

	glBindTexture(GL_TEXTURE_2D, texId[2]);
	glDisable(GL_LIGHTING);
	glColor3d(1, 1, 1);
	glBegin(GL_TRIANGLES);
	for (double u = 0; u <= 0.9001; u += 0.1)
	{
		for (double v = 0; v <= 0.9001; v += 0.1)
		{
			double P1[3] = { 0 };
			double P2[3] = { 0 };
			double P3[3] = { 0 };
			double P4[3] = { 0 };

			P1[0] = f3(4, 4, u, v, fire, 0);
			P1[1] = f3(4, 4, u, v, fire, 1);
			P1[2] = f3(4, 4, u, v, fire, 2);

			P2[0] = f3(4, 4, u + 0.1, v, fire, 0);
			P2[1] = f3(4, 4, u + 0.1, v, fire, 1);
			P2[2] = f3(4, 4, u + 0.1, v, fire, 2);

			P3[0] = f3(4, 4, u, v + 0.1, fire, 0);
			P3[1] = f3(4, 4, u, v + 0.1, fire, 1);
			P3[2] = f3(4, 4, u, v + 0.1, fire, 2);

			P4[0] = f3(4, 4, u + 0.1, v + 0.1, fire, 0);
			P4[1] = f3(4, 4, u + 0.1, v + 0.1, fire, 1);
			P4[2] = f3(4, 4, u + 0.1, v + 0.1, fire, 2);

			glTexCoord2d(u, v);
			glVertex3dv(P1);
			glTexCoord2d(u+0.1, v);
			glVertex3dv(P2);
			glTexCoord2d(u, v+0.1);
			glVertex3dv(P3);

			glTexCoord2d(u+0.1, v);
			glVertex3dv(P2);
			glTexCoord2d(u, v+0.1);
			glVertex3dv(P3);
			glTexCoord2d(u+0.1, v+0.1);
			glVertex3dv(P4);

			if(u == 0.4 && v == 0.4)
			{ 
				if(changefire == 0 && changefire == 0.4 && changefire == 0.8)
				{
					Vector<double> tmp1 = { P1[0], P1[1], P1[2] };
					Vector<double> tmp3 = { P3[0], P3[1], P3[2] };
					track.push_back(tmp1);
					track.push_back(tmp3);
				}
				else
				{
					Vector<double> tmp2 = { P2[0], P2[1], P2[2] };
					Vector<double> tmp4 = { P4[0], P4[1], P4[2] };
					track.push_back(tmp2);
					track.push_back(tmp4);
				}
			}
		}
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void DrowRocket(Point A, Point angle)
{
	double Korp[][3] = {
		{0.5, 0.5, -1.375 },
		{0.5, 0.5, 1.325 },
		{-0.5, 0.5, -1.325 },
		{-0.5, 0.5, 1.325 },
		{-0.5, -0.5, -1.325 },
		{-0.5, -0.5, 1.325 },
		{0.5, -0.5, -1.325 },
		{0.5, -0.5, 1.325 },
		{0, 0, 2.325 }
	};

	for (int i = 0; i < 9; i++)
	{
		Point P = Korp[i];
		P = P.rottX(angle.x);
		P = P.rottZ(angle.z);
		Korp[i][0] = P.x + A.x;
		Korp[i][1] = P.y + A.y;
		Korp[i][2] = P.z + A.z;
	}

	double Engines[][3] = {
		{0.5, -0.25, 0.325},
		{0.5, 0.25, 0.325},
		{0.5, 0.25, -1.625},
		{0.5, -0.25, -1.625},
		{1, -0.25, -1.625},
		{1, 0.25, -1.625},

		{0.25, 0.5, 0.325},
		{-0.25, 0.5, 0.325},
		{-0.25, 0.5, -1.625},
		{0.25, 0.5, -1.625},
		{0.25, 1, -1.625},
		{-0.25, 1, -1.625},

		{-0.5, 0.25, 0.325},
		{-0.5, -0.25, 0.325},
		{-0.5, -0.25, -1.625},
		{-0.5, 0.25, -1.625},
		{-1, 0.25, -1.625},
		{-1, -0.25, -1.625},

		{-0.25, -0.5, 0.325},
		{0.25, -0.5, 0.325},
		{0.25, -0.5, -1.625},
		{-0.25, -0.5, -1.625},
		{-0.25, -1, -1.625},
		{0.25, -1, -1.625},
	};

	DrowFire({ {0.5, 0.25, -1.625}, {0.5, -0.25, -1.625}, {1, 0.25, -1.625}, {1, -0.25, -1.625} }, A, angle);
	DrowFire({ {-0.25, 1, -1.625}, {-0.25, 0.5, -1.625}, {0.25, 1, -1.625}, {0.25, 0.5, -1.625} }, A, angle);
	DrowFire({ {-1, 0.25, -1.625}, {-1, -0.25, -1.625}, {-0.5, 0.25, -1.625}, {-0.5, -0.25, -1.625} }, A, angle);
	DrowFire({ {-0.25, -0.5, -1.625}, {-0.25, -1, -1.625}, {0.25, -0.5, -1.625}, {0.25, -1, -1.625} }, A, angle);

	if (changefire >= 1)
		changefire = 0;
	else changefire += 0.2;

	for (int i = 0; i < 24; i++)
	{
		Point P = Engines[i];
		P = P.rottX(angle.x);
		P = P.rottZ(angle.z);
		Engines[i][0] = P.x + A.x;
		Engines[i][1] = P.y + A.y;
		Engines[i][2] = P.z + A.z;
	}

	glBindTexture(GL_TEXTURE_2D, texId[1]);

	glBegin(GL_QUADS);
	glColor3d(0.5, 0.5, 0.5);
	Point P = Norm(Korp[3], Korp[1], Korp[0]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.1123, 0.01953125);
	glVertex3dv(Korp[0]);
	glTexCoord2d(0.1123, 0.63086);
	glVertex3dv(Korp[1]);
	glTexCoord2d(0.41406, 0.63086);
	glVertex3dv(Korp[3]);
	glTexCoord2d(0.41406, 0.01953125);
	glVertex3dv(Korp[2]);
	glColor3d(0.7, 0.7, 0.7);
	P = Norm(Korp[5], Korp[3], Korp[2]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.1123, 0.01953125);
	glVertex3dv(Korp[2]);
	glTexCoord2d(0.1123, 0.63086);
	glVertex3dv(Korp[3]);
	glTexCoord2d(0.41406, 0.63086);
	glVertex3dv(Korp[5]);
	glTexCoord2d(0.41406, 0.01953125);
	glVertex3dv(Korp[4]);

	P = Norm(Korp[7], Korp[5], Korp[4]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.1123, 0.01953125);
	glVertex3dv(Korp[4]);
	glTexCoord2d(0.1123, 0.63086);
	glVertex3dv(Korp[5]);
	glTexCoord2d(0.41406, 0.63086);
	glVertex3dv(Korp[7]);
	glTexCoord2d(0.41406, 0.01953125);
	glVertex3dv(Korp[6]);

	P = Norm(Korp[1], Korp[7], Korp[6]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.1123, 0.01953125);
	glVertex3dv(Korp[6]);
	glTexCoord2d(0.1123, 0.63086);
	glVertex3dv(Korp[7]);
	glTexCoord2d(0.41406, 0.63086);
	glVertex3dv(Korp[1]);
	glTexCoord2d(0.41406, 0.01953125);
	glVertex3dv(Korp[0]);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(0.7, 0.7, 0.7);
	P = Norm(Korp[0], Korp[6], Korp[4]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.1123, 0.63965);
	glVertex3dv(Korp[0]);
	glTexCoord2d(0.1123, 0.941406);
	glVertex3dv(Korp[6]);
	glTexCoord2d(0.41406, 0.941406);
	glVertex3dv(Korp[4]);
	glTexCoord2d(0.41406, 0.63965);
	glVertex3dv(Korp[2]);
	glEnd();

	for (int i = 0; i < 4; i++)
	{
		glBegin(GL_QUADS);
		glColor3d(0.966, 0.559, 0.31);
		P = Norm(Engines[i * 6 + 4], Engines[i * 6 + 3], Engines[i * 6 + 2]);
		glNormal3d(P.x, P.y, P.z);
		glTexCoord2d(0.63965, 0.58398);
		glVertex3dv(Engines[i * 6 + 2]);
		glTexCoord2d(0.78711, 0.58398);
		glVertex3dv(Engines[i * 6 + 3]);
		glTexCoord2d(0.78711, 0.4375);
		glVertex3dv(Engines[i * 6 + 4]);
		glTexCoord2d(0.63965, 0.4375);
		glVertex3dv(Engines[i * 6 + 5]);

		glColor3d(0.7, 0.0, 0.1);
		P = Norm(Engines[i * 6 + 4], Engines[i * 6 + 5], Engines[i * 6 + 1]);
		glNormal3d(P.x, P.y, P.z);
		glTexCoord2d(0.63965, 0.01953);
		glVertex3dv(Engines[i * 6 + 4]);
		glTexCoord2d(0.78711, 0.01953);
		glVertex3dv(Engines[i * 6 + 5]);
		glTexCoord2d(0.78711, 0.42871);
		glVertex3dv(Engines[i * 6 + 1]);
		glTexCoord2d(0.63965, 0.42871);
		glVertex3dv(Engines[i * 6 + 0]);

		P = Norm(Engines[i * 6 + 0], Engines[i * 6 + 1], Engines[i * 6 + 2]);
		glNormal3d(P.x, P.y, P.z);
		glTexCoord2d(0.63965, 0.42871);
		glVertex3dv(Engines[i * 6 + 0]);
		glTexCoord2d(0.78711, 0.42871);
		glVertex3dv(Engines[i * 6 + 1]);
		glTexCoord2d(0.78711, 0.01953);
		glVertex3dv(Engines[i * 6 + 2]);
		glTexCoord2d(0.63965, 0.01953);
		glVertex3dv(Engines[i * 6 + 3]);
		glEnd();

		glBegin(GL_TRIANGLES);
		P = Norm(Engines[i * 6 + 0], Engines[i * 6 + 3], Engines[i * 6 + 4]);
		glNormal3d(P.x, P.y, P.z);
		glTexCoord2d(0.631836, 0.38574);
		glVertex3dv(Engines[i * 6 + 0]);
		glTexCoord2d(0.631836, 0.01953);
		glVertex3dv(Engines[i * 6 + 3]);
		glTexCoord2d(0.4873, 0.01953);
		glVertex3dv(Engines[i * 6 + 4]);

		P = Norm(Engines[i * 6 + 5], Engines[i * 6 + 2], Engines[i * 6 + 1]);
		glNormal3d(P.x, P.y, P.z);
		glTexCoord2d(0.631836, 0.38574);
		glVertex3dv(Engines[i * 6 + 1]);
		glTexCoord2d(0.631836, 0.01953);
		glVertex3dv(Engines[i * 6 + 2]);
		glTexCoord2d(0.4873, 0.01953);
		glVertex3dv(Engines[i * 6 + 5]);
		glEnd();

		//glRotated(90, 0, 0, 1);
	}

	glBegin(GL_TRIANGLES);
	glColor3d(0.8, 0.1, 0.2);
	P = Norm(Korp[8], Korp[1], Korp[3]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.72852, 0.97168);
	glVertex3dv(Korp[8]);
	glTexCoord2d(0.58, 0.63965);
	glVertex3dv(Korp[1]);
	glTexCoord2d(0.87695, 0.63965);
	glVertex3dv(Korp[3]);

	P = Norm(Korp[8], Korp[3], Korp[5]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.72852, 0.97168);
	glVertex3dv(Korp[8]);
	glTexCoord2d(0.58, 0.63965);
	glVertex3dv(Korp[3]);
	glTexCoord2d(0.87695, 0.63965);
	glVertex3dv(Korp[5]);

	P = Norm(Korp[8], Korp[5], Korp[7]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.72852, 0.97168);
	glVertex3dv(Korp[8]);
	glTexCoord2d(0.58, 0.63965);
	glVertex3dv(Korp[5]);
	glTexCoord2d(0.87695, 0.63965);
	glVertex3dv(Korp[7]);

	P = Norm(Korp[8], Korp[7], Korp[1]);
	glNormal3d(P.x, P.y, P.z);
	glTexCoord2d(0.72852, 0.97168);
	glVertex3dv(Korp[8]);
	glTexCoord2d(0.58, 0.63965);
	glVertex3dv(Korp[7]);
	glTexCoord2d(0.87695, 0.63965);
	glVertex3dv(Korp[1]);
	glEnd();

	if (track.size() >= 200)
	{
		track.erase(track.begin(), track.begin() + 8);
	}

	glColor3d(0.8, 0.8, 0.8);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	if(1)
	{ 
		double max_track = 8;
		for (int i = 0; i < track.size(); i++)
		{
			double P[3] = { track[i][0], track[i][1], track[i][2] };

			glPointSize(max_track + 2 - (double)i / ((double)track.size()/max_track));
			glBegin(GL_POINTS);
			glColor3d(1-(double)i*2/ (track.size() * 3), 1-(double)i*2/ (track.size() * 3), 1-(double)i*2/ (track.size() * 3));
			glVertex3dv(P);
			glEnd();
		}
	}
	else
	{
		glPointSize(6);
		glBegin(GL_POINTS);
		for (int i = 0; i < track.size(); i++)
		{
			double P[3] = { track[i][0], track[i][1], track[i][2] };

			glColor3d(1 - (double)i * 2 / 600, 1 - (double)i * 2 / 600, 1 - (double)i * 2 / 600);
			glVertex3dv(P);
		}
		glEnd();
	}
	glEnable(GL_LIGHTING);
	glEnable(GL_TEXTURE_2D);

}

Point Norm(double * A, double * B, double * C)
{
	/*Point a = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };
	Point b = { C[0] - A[0], C[1] - A[1], C[2] - A[2] };
	
	double x = a.y * b.z - b.y * a.z;
	double y = -a.x * b.z + b.x * a.z;
	double z = a.x * b.y - b.y * a.y;*/
	
	double Ak = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	double Bk = A[2] * (B[0] - C[0]) + B[2] * (C[0] - A[0]) + C[2] * (A[0] - B[0]);
	double Ck = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);

	double x = Ak;
	double y = Bk;
	double z = Ck;

//	if (norm_from_the_cam)
//	{
//		double Dk = A[0] * (B[1] * C[2] - C[1] * B[2]) + B[0] * (C[1] * A[2] - A[1] * C[2]) + C[0] * (A[1] * B[2] - B[1] * A[2]);
//		Dk *= -1;
//		double S = Ak * camera.pos.X() + Bk * camera.pos.Y() + Ck * camera.pos.Z() + Dk;
//		if (S <= 0)
//		{
//			x *= -1;
//			y *= -1;
//			z *= -1;
//		}
//	}

	double l = sqrt(x * x + y * y + z * z);
	x /= l;
	y /= l;
	z /= l;

	return {x,y,z};
}

//double Newtext(double k, int i)
//{
//	if (i == 0)
//		return k / 1610.;
//	else
//		return (1610 - k) / 1610.;
//}



