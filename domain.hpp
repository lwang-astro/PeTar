#pragma once

inline void CalculateBoundaryOfDomain(const PS::S32 &np,
				      const PS::F64vec pos_sample[],
				      const PS::S32 cid,
				      const PS::S32 &istart,
				      const PS::S32 &iend,
				      const PS::F64ort pos_root_domain,
				      PS::F64 & xlow,
				      PS::F64 & xhigh) {
    if(istart == 0) {
	xlow  = pos_root_domain.low_[cid];
    } 
    else {
	xlow  = 0.5 * (pos_sample[istart-1][cid] + pos_sample[istart][cid]);
    }
    if(iend == np - 1) {
	xhigh = pos_root_domain.high_[cid];
    } 
    else {
	xhigh = 0.5 * (pos_sample[iend][cid] + pos_sample[iend+1][cid]);
    }
}

template<class Tsys>
inline void DomainDecision(PS::DomainInfo & dinfo,
			   const Tsys & system){
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    assert(n_proc%4 == 0);
    if(n_proc%4 == 0){
	int n_loc_tmp = system.getNumberOfParticleLocal();
	int * n_recv_tmp = new int[n_proc];
	int * n_recv_disp_tmp = new int[n_proc+1];
	PS::F64vec * pos_loc_tmp = new PS::F64vec[n_loc_tmp];
	for(int i=0; i<n_loc_tmp; i++){
	    pos_loc_tmp[i].z = system[i].pos.z;
	    if(system[i].pos.x * system[i].pos.y > 0.0){
		pos_loc_tmp[i].x = fabs(system[i].pos.x);
		pos_loc_tmp[i].y = fabs(system[i].pos.y);
	    }
	    else{
		pos_loc_tmp[i].x = fabs(system[i].pos.y);
		pos_loc_tmp[i].y = fabs(system[i].pos.x);
	    }
	}
	PS::Comm::allGather(&n_loc_tmp, 1, n_recv_tmp);
	n_recv_disp_tmp[0] = 0;
	for(int i=0; i<n_proc; i++){
	    n_recv_disp_tmp[i+1] = n_recv_disp_tmp[i] + n_recv_tmp[i];
	}
	int n_glb_tmp = n_recv_disp_tmp[n_proc];
	PS::F64vec * pos_glb_tmp = new PS::F64vec[n_glb_tmp];
	PS::Comm::allGatherV(pos_loc_tmp, n_loc_tmp, pos_glb_tmp, n_recv_tmp, n_recv_disp_tmp);
	PS::S32 n_proc_quat = n_proc / 4;
	fout_debug<<"n_proc, n_proc_quat="<<n_proc<<" "<<n_proc_quat<<std::endl;
	PS::S32 nx_quat = sqrt((PS::F64)n_proc_quat-0.000001)+1;
	while( n_proc_quat % nx_quat != 0) nx_quat++;
	PS::S32 ny_quat = n_proc_quat / nx_quat;
	PS::S32 nx = nx_quat*2;
	PS::S32 ny = ny_quat*2;
	PS::S32 nz = 1;
	if(PS::Comm::getRank() == 0){
	    fout_debug<<"nx_quat, ny_quat, nx, ny, nz= "
		      <<nx_quat<<" "<<ny_quat<<" "<<nx<<" "<<ny<<" "<<nz<<std::endl;
	}
	dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
	PS::F64ort * pos_domain_tmp = new PS::F64ort[n_proc];
	if(PS::Comm::getRank() == 0){
	    fout_debug<<"n_glb_tmp="<<n_glb_tmp<<std::endl;
	    PS::S32 * istart = new PS::S32[n_proc_quat];
	    PS::S32 * iend   = new PS::S32[n_proc_quat];
	    PS::F64ort pos_root_domain_tmp = PS::F64ort( PS::F64vec(0.0, 0.0, -PS::LARGE_FLOAT), PS::F64vec(PS::LARGE_FLOAT, PS::LARGE_FLOAT, PS::LARGE_FLOAT));
#ifdef __HPC_ACE__
	    std::sort(pos_glb_tmp, pos_glb_tmp+n_glb_tmp, PS::LessOPX());
#else
	    std::sort(pos_glb_tmp, pos_glb_tmp+n_glb_tmp, 
		      [](const PS::F64vec & l, const PS::F64vec & r)->bool{return l.x < r.x;}
		      );
#endif
	    fout_debug<<"pos_glb_tmp[n_glb_tmp-1].x="<<pos_glb_tmp[n_glb_tmp-1].x<<std::endl;
	    for(PS::S32 i = 0; i < n_proc_quat; i++) {
		istart[i] = ((PS::S64)(i) * (PS::S64)(n_glb_tmp)) / (PS::S64)(n_proc_quat);
		if(i > 0) iend[i-1] = istart[i] - 1;
	    }
	    iend[n_proc_quat-1] = n_glb_tmp - 1;
	    for(PS::S32 ix = 0; ix<nx_quat; ix++) {
		PS::S32 ix0 =  ix      * ny_quat * nz;
		PS::S32 ix1 = (ix + 1) * ny_quat * nz;
		PS::F64 x0, x1;
		CalculateBoundaryOfDomain(n_glb_tmp, pos_glb_tmp, 0, istart[ix0], iend[ix1-1], pos_root_domain_tmp, x0, x1);
		for(PS::S32 i=0; i<ny_quat*nz; i++) {
		    PS::S32 offset = (nx_quat+ix)*ny;
		    pos_domain_tmp[offset+i].low_[0]  = x0;
		    pos_domain_tmp[offset+i].high_[0] = x1;
		    pos_domain_tmp[offset+ny_quat+i].low_[0]  = x0;
		    pos_domain_tmp[offset+ny_quat+i].high_[0] = x1;

		    offset = (nx_quat-1-ix)*ny;
		    pos_domain_tmp[offset+i].low_[0]  = -x1;
		    pos_domain_tmp[offset+i].high_[0] = -x0;
		    pos_domain_tmp[offset+ny_quat+i].low_[0]  = -x1;
		    pos_domain_tmp[offset+ny_quat+i].high_[0] = -x0;
		}
	    }

	    for(PS::S32 ix = 0; ix<nx_quat; ix++) {
		PS::S32 ix0 =  ix      * ny_quat * nz;
		PS::S32 ix1 = (ix + 1) * ny_quat * nz;
#ifdef __HPC_ACE__
	    std::sort(pos_glb_tmp+istart[ix0], pos_glb_tmp+iend[ix1-1]+1, PS::LessOPY());
#else
		std::sort(pos_glb_tmp+istart[ix0], pos_glb_tmp+iend[ix1-1]+1,
			  [](const PS::F64vec & l, const PS::F64vec & r)->bool{return l.y < r.y;}
			  );
#endif
		PS::S32 n_tmp_y = iend[ix1-1] - istart[ix0] + 1;
		for(PS::S32 iy = 0; iy<ny_quat; iy++) {
		    PS::S32 iy0 = ix0 +  iy      * nz;
		    PS::S32 iy1 = ix0 + (iy + 1) * nz;
		    PS::F64 y0, y1;
		    CalculateBoundaryOfDomain(n_tmp_y, pos_glb_tmp+istart[ix0], 1, istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], pos_root_domain_tmp, y0, y1);
		    for(PS::S32 i=0; i<nz; i++) {
			PS::S32 offset = (nx_quat+ix)*ny + ny_quat + iy;
			pos_domain_tmp[offset+i].low_[1]  = y0;
			pos_domain_tmp[offset+i].high_[1] = y1;
			offset = (nx_quat-ix-1)*ny + ny_quat + iy;
			pos_domain_tmp[offset+i].low_[1]  = y0;
			pos_domain_tmp[offset+i].high_[1] = y1;
			offset = (nx_quat+ix)*ny + ny_quat-iy-1;
			pos_domain_tmp[offset+i].low_[1]  = -y1;
			pos_domain_tmp[offset+i].high_[1] = -y0;
			offset = (nx_quat-ix-1)*ny + ny_quat-iy-1;
			pos_domain_tmp[offset+i].low_[1]  = -y1;
			pos_domain_tmp[offset+i].high_[1] = -y0;
		    }
		}
	    }
	    delete [] istart;
	    delete [] iend;
	}
	PS::Comm::broadcast(pos_domain_tmp, n_proc);
	for(PS::S32 i=0; i<n_proc; i++){
	    pos_domain_tmp[i].low_.z  = -PS::LARGE_FLOAT;
	    pos_domain_tmp[i].high_.z =  PS::LARGE_FLOAT;
	    dinfo.setPosDomain(i, pos_domain_tmp[i]);
	}
	delete [] n_recv_tmp;
	delete [] n_recv_disp_tmp;
	delete [] pos_loc_tmp;
	delete [] pos_glb_tmp;
	delete [] pos_domain_tmp;
    }
}
