class KickDriftStep{
    typedef std::array<double,2> KDPair;

    double ds_;  // full step size
    int mode_;   // 0: Kick first; 1: drift first
    int count_one_step_;  // counts for one step case
    int count_continue_;  // counts for continue case
    std::vector<KDPair> coff_one_step_; // cofficient table for one full step
    std::vector<KDPair> coff_continue_; // cofficient table for continuing step (mergin first and last step)

public:

    KickDriftStep(const double _ds): ds_(_ds), mode_(0), count_one_step_(0), count_continue_(0) {
#if (defined KDKDK_2ND) || (defined KDKDK_4TH)
        coff_one_step_.reserve(3);
        coff_one_step_.push_back(KDPair({1.0/6.0, 0.5}));
        coff_one_step_.push_back(KDPair({4.0/6.0, 0.5}));
        coff_one_step_.push_back(KDPair({1.0/6.0, 0.0}));

        coff_continue_.reserve(2);
        coff_continue_.push_back(KDPair({2.0/6.0, 0.5}));
        coff_continue_.push_back(KDPair({4.0/6.0, 0.5}));
#else
        coff_one_step_.reserve(2);
        coff_one_step_.push_back(KDPair({0.5, 1.0}));
        coff_one_step_.push_back(KDPair({0.5, 0.0}));
        assert(coff_one_step_.size()==2);

        coff_continue_.reserve(1);
        coff_continue_.push_back(KDPair({1.0, 1.0}));
        assert(coff_continue_.size()==1);
#endif        
    }

    //! reset step count for one full step case
    void resetCountOneStep() {
        count_one_step_ = 0;
    }

    //! reset step count for continuing case
    void resetCountContinue() {
        count_continue_ = 0;
    }

    //! advance step count for one full step case
    void nextOneStep() {
        count_one_step_++;
        count_one_step_ %= coff_one_step_.size();
    }

    //! advance step count for continuing case
    void nextContinue() {
        count_continue_++;
        count_continue_ %= coff_continue_.size();
    }

    //! get step count for one full step case
    PS::S32 getCountOneStep() const {
        return count_one_step_;
    }

    //! get step count for continue case
    PS::S32 getCountContinue() const {
        return count_continue_;
    }

    //! set DKD mode (default is KDK)
    void setDKDMode() {
        if(count_one_step_||count_continue_) {
            std::cerr<<"Error: in the middel step, cannot switch mode!\n";
            abort();
        }
        mode_ = 1;
    }

    //! set DKD mode (default is KDK)
    void setKDKMode() {
        if(count_one_step_||count_continue_) {
            std::cerr<<"Error: in the middel step, cannot switch mode!\n";
            abort();
        }
        mode_ = 0;
    }

    //! set base step size
    void setStep(const PS::F64 _ds) {
        if(count_one_step_||count_continue_) {
            std::cerr<<"Error: in the middel step, cannot reset step size!\n";
            abort();
        }
        ds_ = _ds;
    }

    //! get base step size
    PS::F64 getStep() const {
        return ds_;
    }
    
    //! Get kick step size for one full step
    PS::F64 getDtKickOneStep() const {
        return coff_one_step_[count_one_step_][mode_]*ds_;
    }

    //! Get drift step size for one full step   
    PS::F64 getDtDriftOneStep() const {
        return coff_one_step_[count_one_step_][1-mode_]*ds_;
    }

    //! Get kick step size for continue case
    PS::F64 getDtKickContinue() const {
        return coff_continue_[count_continue_][mode_]*ds_;
    }

    //! Get drift step size for continue case
    PS::F64 getDtDriftContinue() const {
        return coff_continue_[count_continue_][1-mode_]*ds_;
    }

    //! Get staring step size for continue case
    PS::F64 getDtStartContinue() const {
        if (count_continue_>0) {
            std::cerr<<"Error: try to obtain the first step dt during the middle (count = "<<count_continue_<<")!\n";
            abort();
        }
        return coff_one_step_[0][mode_]*ds_;
    }

    //! Get ending step size for continue case
    PS::F64 getDtEndContinue() const {
        if (count_continue_>0) {
            std::cerr<<"Error: try to obtain the last step dt during the middle (count = "<<count_continue_<<")!\n";
            abort();
        }
        return coff_one_step_[coff_one_step_.size()-1][mode_]*ds_;
    }
};
