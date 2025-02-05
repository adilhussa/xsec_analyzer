#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class CC1mu1cohpiBinScheme : public BinSchemeBase {

  public:

    CC1mu1cohpiBinScheme();
    virtual void DefineBlocks() override;
};
