# server-specific valuables
include LSDMap.inc

# source directory
src_path=src

# RMSD code
rmsd_f90=src/comp_rmsd_3d.f90

lib=${libpath} ${lib_parpack_arpack} ${lib_lapack_blas}
mpiff=${mpif90} ${fflags}
ff=${f90} ${fflags}
c=${cc} ${cflags}

#parallel:
#	${mpiff} ${src_path}p_wlsdmap_embed.f90 ${rmsd_f90} -o p_wlsdmap_embed ${lib}

default: parallel serial

parallel: main split_rmsd p_wlsdmap_embed
	mkdir -p input err embed
serial: s_rmsd_neighbor s_wlsdmap
	mkdir -p rmsd err

main: ${src_path}/data.o ${src_path}/parallel.o ${src_path}/p_wlsdmap.o ${src_path}/qsort.o ${src_path}/ftn_c.o ${src_path}/p_local_mds.o ${src_path}/qcprot.o  ${src_path}/p_rmsd_neighbor.o  ${src_path}/main.o
	${mpiff} $^ -o $@ ${lib}
split_rmsd: ${src_path}/split_rmsd.f90
	${ff} $^ -o $@
p_wlsdmap_embed:  ${src_path}/ftn_c.o ${src_path}/p_wlsdmap_embed.o ${src_path}/qcprot.o
	${mpiff} $^ -o $@ ${lib}
s_rmsd_neighbor: ${src_path}/qsort.o ${src_path}/ftn_c.o ${src_path}/s_rmsd_neighbor.o ${src_path}/qcprot.o
	${ff} $^ -o $@
s_wlsdmap: ${src_path}/s_wlsdmap.o
	${mpiff} $^ -o $@ ${lib}


${src_path}/main.o: ${src_path}/main.f90
	${mpiff} -c $^ -o $@
${src_path}/p_wlsdmap_embed.o: ${src_path}/p_wlsdmap_embed.f90
	${mpiff} -c $^ -o $@
${src_path}/p_wlsdmap.o: ${src_path}/p_wlsdmap.f90
	${mpiff} -c $^ -o $@
${src_path}/p_local_mds.o : ${src_path}/p_local_mds.f90
	${mpiff} -c $^ -o $@
${src_path}/p_rmsd_neighbor.o : ${src_path}/p_rmsd_neighbor.f90
	${mpiff} -c $^ -o $@
${src_path}/s_wlsdmap.o: ${src_path}/s_wlsdmap.f90
	${ff} -c $^ -o $@
${src_path}/s_rmsd_neighbor.o: ${src_path}/s_rmsd_neighbor.f90
	${ff} -c $^ -o $@
${src_path}/qsort.o: ${src_path}/qsort.f90
#	cd ${src_path}; ${ff} -c qsort.f90; cd ../
	${ff} -c $^ -o $@
${src_path}/ftn_c.o: ${src_path}/ftn_c.f90
	${ff} -c $^ -o $@
${src_path}/qcprot.o: ${src_path}/qcprot.c
	${c} -c $^ -o $@
${src_path}/data.o: ${src_path}/data.f90
	${mpiff} -c $^ -o $@
${src_path}/parallel.o: ${src_path}/parallel.f90
	${mpiff} -c $^ -o $@

clean:
	rm -f p_rmsd_neighbor p_local_mds p_wlsdmap p_wlsdmap_embed split_rmsd s_rmsd_neighbor s_wlsdmap *.mod src/*.o
