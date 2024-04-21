#pragma once
#ifdef BSE_BASE
#include "bse_interface.h"
#endif
#ifdef NAN_CHECK_DEBUG
#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif
#endif

#ifdef BSE_BASE
enum class BinaryInterruptState:int {none = 0, form = 1, type_change = 2, start_roche = 3, end_roche = 4, contact = 5, start_symbiotic = 6, end_symbiotic = 7, common_envelope = 8 , giant = 9, collision = 10, blue_straggler = 11, no_remain = 12, disrupt = 13, tide = 14};
#else
enum class BinaryInterruptState:int {none = 0, form = 1, exchange = 2, collision = 3, delaycollision = 4};
#endif
#define BINARY_STATE_ID_SHIFT 4
#define BINARY_INTERRUPT_STATE_MASKER 0xF

/// Basic particle class
class ParticleBase{
public:
    // necessary variables, should not be touched
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::S64 binary_state; // contain two parts, low bits (first BINARY_STATE_ID_SHIFT bits) is binary interrupt state and high bits are pair ID
#ifdef STELLAR_EVOLUTION
    // for stellar evolution
    PS::F64 radius;
    PS::F64 dm;
    PS::F64 time_record; 
    PS::F64 time_interrupt;
#ifdef BSE_BASE
    StarParameter star; // SSE/BSE based package stellar parameters
#endif
#endif

    //! save pair id in binary_state with shift bit size of BINARY_STATE_ID_SHIFT
    void setBinaryPairID(const PS::S64 _id) {
        binary_state = (binary_state&BINARY_INTERRUPT_STATE_MASKER) | (_id<<BINARY_STATE_ID_SHIFT);
    }

    //! save binary interrupt state in the first  BINARY_STATE_ID_SHIFT bit in binary_state
    void setBinaryInterruptState(const BinaryInterruptState _state) {
        binary_state = ((binary_state>>BINARY_STATE_ID_SHIFT)<<BINARY_STATE_ID_SHIFT) | int(_state);
    }

    //! get binary interrupt state from binary_state
    BinaryInterruptState getBinaryInterruptState() const {
        return static_cast<BinaryInterruptState>(binary_state&BINARY_INTERRUPT_STATE_MASKER);
    }

    //! get pair ID from binary_state 
    PS::S64 getBinaryPairID() const {
        return (binary_state>>BINARY_STATE_ID_SHIFT);
    }

    // -------------------------

    //! defaulted constructor 
    ParticleBase(): mass(0.0) { 
        binary_state = 0;
#ifdef STELLAR_EVOLUTION
        radius = 0.0;
        dm = 0.0;
        time_record = 0.0;
        time_interrupt = 0.0;
#endif
    }

    template<class Tp>
    ParticleBase(const Tp &p) { DataCopy(p);}

    //! constructor 
    ParticleBase(const PS::F64 _mass, 
                 const PS::F64vec & _pos, 
                 const PS::F64vec & _vel): mass(_mass), pos(_pos), vel(_vel) {
        binary_state = 0;
#ifdef STELLAR_EVOLUTION
        radius = 0.0;
        dm = 0.0;
        time_record = 0.0;
        time_interrupt = 0.0;
#ifdef BSE_BASE
        star.initial(0.0);
#endif
#endif
    }

#ifdef STELLAR_EVOLUTION
#ifdef BSE_BASE
    //! constructor 
    ParticleBase(const PS::F64 _mass, const PS::F64vec & _pos, const PS::F64vec & _vel, const PS::S64 _binary_state,
                 const PS::F64 _radius, const PS::F64 _dm, 
                 const PS::F64 _time_record, const PS::F64 _time_interrupt, const StarParameter& _star): 
        mass(_mass), pos(_pos), vel(_vel), binary_state(_binary_state), radius(_radius), dm(_dm), 
        time_record(_time_record), time_interrupt(_time_interrupt), star(_star) {}
#else
    //! constructor 
    ParticleBase(const PS::F64 _mass, const PS::F64vec & _pos, const PS::F64vec & _vel, const PS::S64 _binary_state,
                 const PS::F64 _radius, const PS::F64 _dm, 
                 const PS::F64 _time_record, const PS::F64 _time_interrupt):
        mass(_mass), pos(_pos), vel(_vel), binary_state(_binary_state), radius(_radius), dm(_dm), 
        time_record(_time_record), time_interrupt(_time_interrupt) {}
#endif
#else
    //! constructor 
    ParticleBase(const PS::F64 _mass, const PS::F64vec & _pos, const PS::F64vec & _vel, const PS::S64 _binary_state):
        mass(_mass), pos(_pos), vel(_vel), binary_state(_binary_state) {}
#endif

#ifdef STELLAR_EVOLUTION
    //! write class data with ASCII format
    /*! @param[in] _fout: file IO for write
     */
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %lld %26.17e %26.17e %26.17e %26.17e ",
                this->mass, 
                this->pos.x, this->pos.y, this->pos.z,  
                this->vel.x, this->vel.y, this->vel.z,
                this->binary_state, this->radius, this->dm, this->time_record, this->time_interrupt);
#ifdef BSE_BASE
        star.writeAscii(fp);
#endif
    }

    //! read class data with ASCII format
    /*! @param[in] _fin: file IO for read
     */
    void readAscii(FILE* fp) {
        PS::S64 rcount=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lld %lf %lf %lf %lf",
                              &this->mass, 
                              &this->pos.x, &this->pos.y, &this->pos.z,
                              &this->vel.x, &this->vel.y, &this->vel.z,
                              &this->binary_state,
                              &this->radius, &this->dm, &this->time_record, &this->time_interrupt);
        if(rcount<12) {
            std::cerr<<"Error: Data reading fails! requiring data number is 12, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
#ifdef BSE_BASE
        star.readAscii(fp);
#endif
    }
#else
    //! write class data with ASCII format
    /*! @param[in] _fout: file IO for write
     */
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %lld",
                this->mass, 
                this->pos.x, this->pos.y, this->pos.z,  
                this->vel.x, this->vel.y, this->vel.z,
                this->binary_state);
    }

    //! read class data with ASCII format
    /*! @param[in] _fin: file IO for read
     */
    void readAscii(FILE* fp) {
        PS::S64 rcount=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lld",
                              &this->mass, 
                              &this->pos.x, &this->pos.y, &this->pos.z,
                              &this->vel.x, &this->vel.y, &this->vel.z,
                              &this->binary_state);
        if(rcount<8) {
            std::cerr<<"Error: Data reading fails! requiring data number is 7, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
    }
#endif

    //! write class data with BINARY format
    /*! @param[in] _fout: file IO for write
     */
    void writeBinary(FILE* fp) const{
        fwrite(&(this->mass), sizeof(ParticleBase), 1, fp);
    }


    //! read class data with BINARY format
    /*! @param[in] _fin: file IO for read
     */
    void readBinary(FILE* fp) {
        size_t rcount=fread(&(this->mass), sizeof(ParticleBase), 1, fp);
        if(rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            std::cerr<<"Check your input data, whether the consistent features (interrupt mode and external mode) are used in configuring petar and the data generation\n";
            abort();
        }
    }

    //! for print debugging
    void print(std::ostream & fout) const{
        fout<<" mass="<<mass
            <<" pos="<<pos
            <<" vel="<<vel
            <<" binary_state="<<binary_state;
#ifdef STELLAR_EVOLUTION
        fout<<" radius="<<radius
            <<" dm="<<dm
            <<" time_record="<<time_record
            <<" time_interrupt="<<time_interrupt;
#ifdef BSE_BASE
        star.print(fout);
#endif
#endif
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"mass"
             <<std::setw(_width)<<"pos.x"
             <<std::setw(_width)<<"pos.y"
             <<std::setw(_width)<<"pos.z"
             <<std::setw(_width)<<"vel.x"
             <<std::setw(_width)<<"vel.y"
             <<std::setw(_width)<<"vel.z"
             <<std::setw(_width)<<"bin_stat";
#ifdef STELLAR_EVOLUTION
        _fout<<std::setw(_width)<<"radius"
             <<std::setw(_width)<<"dm"
             <<std::setw(_width)<<"t_record"
             <<std::setw(_width)<<"t_interrupt";
#ifdef BSE_BASE
        StarParameter::printColumnTitle(_fout, _width);
#endif
#endif
    }

    //! print column title with meaning (each line for one column)
    /*! @param[out] _fout: std::ostream output object
      @param[in] _counter: offset of the number counter for each line to indicate the column index (defaulted 0)
      @param[in] _offset: the printing whitespace offset for each line (defaulted 0)
      \return: the total counter of columns
     */
    static int printTitleWithMeaning(std::ostream & _fout, const int _counter=0, const int _offset=0) {
        int counter = _counter;
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". mass: mass of particle\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<"-"<<counter+2<<". pos.[x/y/z]: 3D position of particle\n";
        counter+=3;
        _fout<<std::setw(_offset)<<" "<<counter<<"-"<<counter+2<<". vel.[x/y/z]: 3D velocity of particle\n";
        counter+=3;
        _fout<<std::setw(_offset)<<" "<<counter<<". bin_stat: binary status storing pair id and status [formatted] (0.0)\n";
#ifdef STELLAR_EVOLUTION
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". radius: stellar radius for merger checker (0.0)\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". dm: mass change (0.0)\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". t_record: time record of last check (0.0)\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". t_interrupt: time for next evolution check (0.0)\n";
#ifdef BSE_BASE
        counter = StarParameter::printTitleWithMeaning(_fout, counter, _offset);
#endif
#endif
        return counter;
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20) const{
        _fout<<std::setw(_width)<<mass
             <<std::setw(_width)<<pos.x
             <<std::setw(_width)<<pos.y
             <<std::setw(_width)<<pos.z
             <<std::setw(_width)<<vel.x
             <<std::setw(_width)<<vel.y
             <<std::setw(_width)<<vel.z
             <<std::setw(_width)<<binary_state;
#ifdef STELLAR_EVOLUTION
        _fout<<std::setw(_width)<<radius
             <<std::setw(_width)<<dm
             <<std::setw(_width)<<time_record
             <<std::setw(_width)<<time_interrupt;
#ifdef BSE_BASE
        star.printColumn(_fout, _width);
#endif
#endif
    }

    //! print data of class members with pos and vel offset using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[in] _pcm: particle data with position and velocity offset that are added when print data
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
        template <class Tpcm>
        void printColumnWithOffset(Tpcm& _pcm, std::ostream & _fout, const int _width=20) const{
        _fout<<std::setw(_width)<<mass
             <<std::setw(_width)<<pos.x + _pcm.pos.x
             <<std::setw(_width)<<pos.y + _pcm.pos.y
             <<std::setw(_width)<<pos.z + _pcm.pos.z
             <<std::setw(_width)<<vel.x + _pcm.vel.x
             <<std::setw(_width)<<vel.y + _pcm.vel.y
             <<std::setw(_width)<<vel.z + _pcm.vel.z
             <<std::setw(_width)<<binary_state;
#ifdef STELLAR_EVOLUTION
        _fout<<std::setw(_width)<<radius
             <<std::setw(_width)<<dm
             <<std::setw(_width)<<time_record
             <<std::setw(_width)<<time_interrupt;
#ifdef BSE_BASE
        star.printColumn(_fout, _width);
#endif
#endif
    }
    

    //! Copy from another ParticleBase 
    /*! This is used for data transfer between nodes and between soft and hard parts
      @param[in] din: data need to be copied
     */
    template<class Tp>
    void DataCopy(const Tp & din) {
        mass = din.mass;
        pos  = din.pos;
        vel  = din.vel;
        binary_state  = din.binary_state;
#ifdef STELLAR_EVOLUTION
        radius= din.radius;
        dm  = din.dm;
        time_record  = din.time_record;
        time_interrupt = din.time_interrupt;
#ifdef BSE_BASE
        star = din.star;
#endif
#endif
    }

    //! Get mass (required for \ref ARC::chain)
    /*! \return mass
     */
    PS::F64 getMass() {
        return mass;
    }
  
    //! Get position (required for \ref ARC::chain)
    /*! \return position vector (PS::F64[3])
     */
    PS::F64* getPos() {
        return &pos.x;
    }

    //! Get velocity (required for \ref ARC::chain)
    /*! \return velocity vector (PS::F64[3])
     */
    PS::F64* getVel() {
        return &vel.x;
    }

    //!Set position (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] x: particle position in x axis
      @param [in] y: particle position in y axis
      @param [in] z: particle position in z axis
    */
    void setPos(const PS::F64 x, const PS::F64 y, const PS::F64 z) {
#ifdef NAN_CHECK_DEBUG
        NAN_CHECK(x);
        NAN_CHECK(y);
        NAN_CHECK(z);
#endif    
        pos.x = x;
        pos.y = y;
        pos.z = z;
    }

    //!Set position (used in soft part)
    void setPos(const PS::F64vec & _pos) {pos = _pos;}
    
    //!Set velocity (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] vx: particle velocity in x axis
      @param [in] vy: particle velocity in y axis 
      @param [in] vz: particle velocity in z axis 
    */
    void setVel(const PS::F64 vx, const PS::F64 vy, const PS::F64 vz) {
#ifdef NAN_CHECK_DEBUG
        NAN_CHECK(vx);
        NAN_CHECK(vy);
        NAN_CHECK(vz);
#endif    
    
        vel.x = vx;
        vel.y = vy;
        vel.z = vz;
    }

    //!Set velocity
    void setVel(const PS::F64vec & _vel) {vel = _vel;};

    //!Set mass (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] m: particle mass
    */
    void setMass(const PS::F64 m) {
#ifdef NAN_CHECK_DEBUG
        NAN_CHECK(m);
#endif

        mass = m;
    }
};
