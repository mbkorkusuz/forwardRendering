#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

class Triangle
{
public:
    int vertexIds[3];

    Triangle();
    Triangle(int vid1, int vid2, int vid3);
    Triangle(const Triangle &other);
};

#endif