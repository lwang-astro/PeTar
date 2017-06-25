#pragma once

/// Basic particle class
class ParticleBase{
public:
    // necessary variables, should not be touched
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    // -------------------------

    //! defaulted constructor 
    ParticleBase() {}

    template<class Tp>
    ParticleBase(const Tp &p) { DataCopy(p);}

    //! constructor 
    ParticleBase(const PS::F64 _mass, 
                 const PS::F64vec & _pos, 
                 const PS::F64vec & _vel): mass(_mass), pos(_pos), vel(_vel) {}

    //! for data write
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
                this->mass, 
                this->pos.x, this->pos.y, this->pos.z,  
                this->vel.x, this->vel.y, this->vel.z);
    }

    //! for data read
    void readAscii(FILE* fp) {
        PS::S64 rcount=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf ",
                              &this->mass, 
                              &this->pos.x, &this->pos.y, &this->pos.z,
                              &this->vel.x, &this->vel.y, &this->vel.z);
        if(rcount<7) {
            std::cerr<<"Error: Data reading fails! requiring data number is 7, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! for dump debugging
    void dump(std::ofstream & fout){
        fout<<" mass="<<mass
            <<" pos="<<pos
            <<" vel="<<vel;
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
    }

    //! Get mass (required for \ref ARC::chain)
    /*! \return mass
     */
    const PS::F64 getMass() const{
        return mass;
    }
  
    //! Get position (required for \ref ARC::chain)
    /*! \return position vector (PS::F64[3])
     */
    const PS::F64* getPos() const{
        return &pos[0];
    }

    //! Get velocity (required for \ref ARC::chain)
    /*! \return velocity vector (PS::F64[3])
     */
    const PS::F64* getVel() const{
        return &vel[0];
    }

    //!Set position (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] x: particle position in x axis
      @param [in] y: particle position in y axis
      @param [in] z: particle position in z axis
    */
    void setPos(const PS::F64 x, const PS::F64 y, const PS::F64 z) {
        NAN_CHECK(x);
        NAN_CHECK(y);
        NAN_CHECK(z);
    
        pos[0] = x;
        pos[1] = y;
        pos[2] = z;
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
        NAN_CHECK(vx);
        NAN_CHECK(vy);
        NAN_CHECK(vz);
    
        vel[0] = vx;
        vel[1] = vy;
        vel[2] = vz;
    }

    //!Set velocity
    void setVel(const PS::F64vec & _vel) {vel = _vel;};

    //!Set mass (required for \ref ARC::chain)
    /*! NAN check will be done
      @param [in] m: particle mass
    */
    void setMass(const PS::F64 m) {
        NAN_CHECK(m);

        mass = m;
    }
};
