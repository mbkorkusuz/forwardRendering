#ifndef __TRANSLATION_H__
#define __TRANSLATION_H__

class Translation
{
public:
    int translationId;
    double tx, ty, tz;

    Translation();
    Translation(int translationId, double tx, double ty, double tz);
};

#endif
