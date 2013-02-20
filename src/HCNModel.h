#ifndef HCNMODEL_H
#define HCNMODEL_H

class HCNModel {
public:
    HCNModel();

    void generateHCNs() { }
};

class HCNGrowthModel : public HCNModel {
public:
    double startTime, endTime;
    int timePointCount;
};

#endif // HCNMODEL_H
