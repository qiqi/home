#ifndef OPERATOR_TESTADJ_H
#define OPERATOR_TESTADJ_H

#include "field.h"

void Randomize ( Field& u );
void Randomize ( EdgeField& U );

int  TestAdj_main (void);

void TestFaceNormalAdjoint       (AMesh* mesh);
void TestGlobalDivergenceAdjoint (AMesh* mesh);
void TestFaceGradientAdjoint     (AMesh* mesh);
void TestNodeGradientAdjoint     (AMesh* mesh);
void TestDivergenceAdjoint       (AMesh* mesh);
void TestLaplaceAdjoint          (AMesh* mesh);
void TestConvectionAdjoint       (AMesh* mesh);

#endif
