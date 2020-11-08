#ifndef SEPSETPRODUCER_HPP_
#define SEPSETPRODUCER_HPP_


/**
 * @author Joseph Ramsey
 */

class SepsetProducer {

public:

    virtual std::vector<Variable*>  getSepset(Variable* a, Variable* b) = 0;

    virtual bool isCollider(Variable* i, Variable* j, Variable* k) = 0;

    virtual bool isNoncollider(Variable* a, Variable* b, Variable* c) = 0;

    virtual bool isIndependent(Variable* a, Variable* b, std::vector<Variable*>  c) = 0;

    virtual double getPValue() = 0;

    virtual double getScore() = 0;

    virtual std::vector<Variable*> getVariables() = 0;

    virtual void setVerbose(bool verbose) = 0;
};

#endif /* SEPSETPRODUCER_HPP_ */
