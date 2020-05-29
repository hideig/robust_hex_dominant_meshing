#include "TetraMeshTool/Toperator.h"
#include <iterator>

namespace VolumeMesh
{
	/************************************************************************/
	/* Tetrahedral Mesh Topological Operation                               */
	/* Topological Operation                                                */
	/************************************************************************/

	//********************************************FLIP OPERATION*********************************************//

	/** topology operation : flip1-3.
	*	\param hf_: index of the HalfFace which contains the insertPoint
	*   \param insertPoint: the coordinate of the insert point which should be inside the half face hf_
	*   \param tetraVec_(optional): return new tetrahedrons 
	*   \return 0 success.
	*           1 fail.
	*/
	int TMeshToperator::flip13(HalfFaceHandle &hf_, Point &insertPoint, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::runtime_error("No Mesh!");

			if (hf_.idx() < 0 || hf_.idx() >= (int)tmesh->size_halfface())
				throw std::range_error("FLIP1-3: The index of HalfFace handle is out of range!");

			// get the tetrahedron which contains the half face hf_
			TetraHandle tetra = tmesh->handle_to_entity(hf_).hedron_handle();
			// get the points in half face hf_
			TetraMesh::HalfFaceVertexIter hfv_iter;
			PointHandle fpoint[3];
			hfv_iter = tmesh->half_face_vertex_iter(hf_);
			fpoint[0] = tmesh->point_handle(hfv_iter.handle()); ++hfv_iter;
			fpoint[1] = tmesh->point_handle(hfv_iter.handle()); ++hfv_iter;
			fpoint[2] = tmesh->point_handle(hfv_iter.handle());
			// get the point not in half face hf_
			PointHandle top;
			TetraMesh::FaceHalfedgeIter fhe_iter = tmesh->face_half_edge_iter(hf_);
			top = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(fhe_iter.handle()))));

			// delete the original tetrahedron
			if (tmesh->erase_tetrahedron(tetra) == TetraHandle(-1))
				throw std::runtime_error("FLIP1-3: Tetrahedron erase error!");

			// add the three new tetrahedrons
			TetraHandle newTetra[3];
			newTetra[0] = tmesh->add_tetrahedron(tmesh->point(top), tmesh->point(fpoint[0]), tmesh->point(fpoint[1]), insertPoint);
			if (newTetra[0] == TetraHandle(-1))
				throw std::runtime_error("FLIP1-3: Tetrahedron add fail!");
			newTetra[1] = tmesh->add_tetrahedron(tmesh->point(top), tmesh->point(fpoint[0]), insertPoint, tmesh->point(fpoint[2]));
			if (newTetra[1] == TetraHandle(-1))
				throw std::runtime_error("FLIP1-3: Tetrahedron add fail!");
			newTetra[2] = tmesh->add_tetrahedron(tmesh->point(top), insertPoint, tmesh->point(fpoint[1]), tmesh->point(fpoint[2]));
			if (newTetra[2] == TetraHandle(-1))
				throw std::runtime_error("FLIP1-3: Tetrahedron add fail!");

			// save the recover data
			rec_flip13_p[0] = top;
			rec_flip13_p[1] = fpoint[0];
			rec_flip13_p[2] = fpoint[1];
			rec_flip13_p[3] = fpoint[2];
			rec_flip13_old_tetra = tetra;
			for (int i = 0; i < 3; ++i)
				rec_flip13_tetra[i] = newTetra[i];

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra[0]);
				tetraVec_->push_back(newTetra[1]);
				tetraVec_->push_back(newTetra[2]);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}


	/** topology operation : flip3-1.
	*   \param midPoint: the coordinate of the shared point by three boundary faces
	*   \param tetraVec_(optional): return new tetrahedrons 
	*   \return 0 success.
	*           1 fail.
	*/
	int TMeshToperator::flip31(HalfFaceHandle hf[3], PointHandle &midPointHandle, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::runtime_error("No Mesh!");

			/* check if the handles are valid */
			if (midPointHandle.idx() < 0 || midPointHandle.idx() >= (int)tmesh->size_point())
				throw std::range_error("FLIP3-1: The index of Point handle is out of range!");
			if (hf[0].idx() < 0 || hf[0].idx() >= (int)tmesh->size_halfface())
				throw std::range_error("FLIP3-1: The index of HalfFace handle is out of range!");
			if (hf[1].idx() < 0 || hf[1].idx() >= (int)tmesh->size_halfface())
				throw std::range_error("FLIP3-1: The index of HalfFace handle is out of range!");
			if (hf[2].idx() < 0 || hf[2].idx() >= (int)tmesh->size_halfface())
				throw std::range_error("FLIP3-1: The index of HalfFace handle is out of range!");

			// get the tetrahedron which contains the midpoint
			TetraHandle tetras[3];
			tetras[0] = tmesh->handle_to_entity(hf[0]).hedron_handle();
			tetras[1] = tmesh->handle_to_entity(hf[1]).hedron_handle();
			tetras[2] = tmesh->handle_to_entity(hf[2]).hedron_handle();

			// get related points
			PointHandle top;
			PointHandle p[3];
			Point midpoint = tmesh->point(midPointHandle);
			TetraMesh::FaceHalfedgeIter fe_it;
			fe_it = tmesh->face_half_edge_iter(hf[0]);
			top = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(fe_it.handle()))));
			p[0] = tmesh->point_handle(tmesh->from_vertex_handle(fe_it.handle()));
			p[1] = tmesh->point_handle(tmesh->to_vertex_handle(fe_it.handle()));
			fe_it = tmesh->face_half_edge_iter(hf[1]);
			for (; fe_it; ++fe_it)
			{
				p[2] = tmesh->point_handle(tmesh->from_vertex_handle(fe_it.handle()));
				if (p[2] != p[0] && p[2] != p[1])
					break;
			}
			
			// delete the original tetrahedron
			if (tmesh->erase_tetrahedron(tetras[0]) == TetraHandle(-1))
				throw std::runtime_error("FLIP3-1: Tetrahedron erase error!");
			if (tmesh->erase_tetrahedron(tetras[1]) == TetraHandle(-1))
				throw std::runtime_error("FLIP3-1: Tetrahedron erase error!");
			if (tmesh->erase_tetrahedron(tetras[2]) == TetraHandle(-1))
				throw std::runtime_error("FLIP3-1: Tetrahedron erase error!");

			// add s new tetrahedron
			TetraHandle newTetra;
			newTetra = tmesh->add_tetrahedron(tmesh->point(top), tmesh->point(p[0]), tmesh->point(p[1]), tmesh->point(p[2]));
			if (newTetra == TetraHandle(-1))
				throw std::runtime_error("FLIP3-1: Tetrahedron add fail!");

			// save the recover data
			rec_flip31_top = top;
			rec_flip31_face_p[0] = p[0];
			rec_flip31_face_p[1] = p[1];
			rec_flip31_face_p[2] = p[2];
			rec_flip31_old_tetra[0] = tetras[0];
			rec_flip31_old_tetra[1] = tetras[1];
			rec_flip31_old_tetra[2] = tetras[2];
			rec_flip31_tetra = newTetra;


			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	/** topology operation : flip2-3.
	*	\param hf_ index of a HalfFace of the shared face
	*   \param tetraVec_(optional) return new tetrahedrons 
	*   \return 1 invalid.
	*			0 success.
	*/
	int TMeshToperator::flip23(HalfFaceHandle &hf_, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (hf_.idx() < 0 || hf_.idx() >= (int)tmesh->size_halfface())
				throw std::range_error("FLIP2-3: The index of HalfFace handle is out of range!");

			// get the two original tetrahedrons
			TetraHandle tetras[2];
			tetras[0] = TetraHandle(hf_.idx()>>2);
			tetras[1] = TetraHandle(tmesh->handle_to_entity(hf_).opposite_face_handle().idx()>>2);

			if ((tetras[0] < 0 || tetras[0] >= (int)tmesh->size_tetrahedron())||(tetras[1] < 0 || tetras[1] >= (int)tmesh->size_tetrahedron()))
				throw std::range_error("FLIP2-3: The index of Tetrahedron handle is out of range!");

			// get related five points
			Point top, bot, p[3];
			PointHandle toph, both, ph[3];
			VolumeMesh::TetraMesh::FaceHalfedgeIter fe_iter;
			fe_iter = tmesh->face_half_edge_iter(hf_);
			top = tmesh->point(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(fe_iter.handle()))));
			bot = tmesh->point(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(fe_iter.handle())))));
			p[0] = tmesh->point(tmesh->from_vertex_handle(fe_iter.handle()));++fe_iter;
			p[1] = tmesh->point(tmesh->from_vertex_handle(fe_iter.handle()));++fe_iter;
			p[2] =tmesh-> point(tmesh->from_vertex_handle(fe_iter.handle()));

			fe_iter = tmesh->face_half_edge_iter(hf_);
			toph = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(fe_iter.handle()))));
			both = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(fe_iter.handle())))));
			ph[0] = tmesh->point_handle(tmesh->from_vertex_handle(fe_iter.handle()));++fe_iter;
			ph[1] = tmesh->point_handle(tmesh->from_vertex_handle(fe_iter.handle()));++fe_iter;
			ph[2] =tmesh-> point_handle(tmesh->from_vertex_handle(fe_iter.handle()));

			// do flip
			// delete two original tetrahedrons
			if (tmesh->erase_tetrahedron(tetras[0])==TetraHandle(-1))
				throw std::runtime_error("FLIP2-3: Tetrahedron erase fail!");
			if (tmesh->erase_tetrahedron(tetras[1])==TetraHandle(-1))
				throw std::runtime_error("FLIP2-3: Tetrahedron erase fail!");

			// add three new tetrahedrons
			TetraHandle newTetra[3];
			newTetra[0] = tmesh->add_tetrahedron(top, p[0], p[1], bot, &tetras[0]);
			if (newTetra[0] == TetraHandle(-1))
				throw std::runtime_error("FLIP2-3: Tetrahedron add fail!");
			newTetra[1] = tmesh->add_tetrahedron(top, p[0], bot, p[2], &tetras[1]);
			if (newTetra[1] == TetraHandle(-1))
				throw std::runtime_error("FLIP2-3: Tetrahedron add fail!");
			newTetra[2] = tmesh->add_tetrahedron(top, bot, p[1], p[2]);
			if (newTetra[2] == TetraHandle(-1))
				throw std::runtime_error("FLIP2-3: Tetrahedron add fail!");

			// save recover data
			rec_flip23_top = toph;
			rec_flip23_bot = both;
			rec_flip23_p[0] = ph[0];
			rec_flip23_p[1] = ph[1];
			rec_flip23_p[2] = ph[2];
			rec_flip23_old_tetra[0] = tetras[0];
			rec_flip23_old_tetra[1] = tetras[1];
			for (int i = 0; i < 3; ++i)
				rec_flip23_tetra[i] = newTetra[i];
			// return the three new tetrahedrons
			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra[0]);
				tetraVec_->push_back(newTetra[1]);
				tetraVec_->push_back(newTetra[2]);
			}
		}
		catch (std::range_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}


	/** topology operation : flip3-2.
	*	\param he_ index of one of the shared halfEdge
	*   \param tetraVec_(optional) return the new tetrahedrons
	*   \return  0 success.
	*			 1 invalid.
	*/
	int TMeshToperator::flip32(HalfEdgeHandle &he_, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (he_.idx() < 0 || he_.idx() >= (int)tmesh->size_halfedge())
				throw std::range_error("FLIP3-2: The index of HalfFace handle is out of range!");

			// get tetrahedrons which contain the edge
			std::vector<TetraHandle> tetras;
			tmesh->edge_star(he_, tetras);
			if (tetras.size() != 3)
				throw std::runtime_error("FLIP3-2: Edge star size error!");

			// get the five related points
			PointHandle toph, both, ph[3];
			Point top, bot, p[3];
			toph = tmesh->point_handle(tmesh->from_vertex_handle(he_));
			both = tmesh->point_handle(tmesh->to_vertex_handle(he_));
			ph[0] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(he_)));
			ph[2] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(he_))));
			ph[1] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(he_)))));

			top = tmesh->point(toph);
			bot = tmesh->point(both);
			p[0] = tmesh->point(ph[0]);
			p[1] = tmesh->point(ph[1]);
			p[2] = tmesh->point(ph[2]);

			// do flip
			// delete the three old tetrahedrons
			if (tmesh->erase_tetrahedron(tetras[0])==TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron erase fail!");
			if (tmesh->erase_tetrahedron(tetras[1])==TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron erase fail!");
			if (tmesh->erase_tetrahedron(tetras[2])==TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron erase fail!");

			// add the two new tetrahedrons
			TetraHandle newTetra[2];
			newTetra[0] = tmesh->add_tetrahedron(top, p[0], p[1], p[2]);
			if (newTetra[0] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");
			newTetra[1] = tmesh->add_tetrahedron(bot, p[2], p[1], p[0]);
			if (newTetra[1] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");

			// save recover data
			rec_flip32_top = toph;
			rec_flip32_bot = both;
			rec_flip32_p[0] = ph[0];
			rec_flip32_p[1] = ph[1];
			rec_flip32_p[2] = ph[2];
			rec_flip32_old_tetra[0] = tetras[0];
			rec_flip32_old_tetra[1] = tetras[1];
			rec_flip32_old_tetra[2] = tetras[2];
			rec_flip32_tetra[0] = newTetra[0];
			rec_flip32_tetra[1] = newTetra[1];

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra[0]);
				tetraVec_->push_back(newTetra[1]);
			}
		}
		catch (std::range_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	/** topology operation : flip2-2.
	*	\param he_ index of the shared halfedge of two coplanar halffaces
	*   \param tetraVec_(optional) return the new tetrahedrons
	*   \return  0 success.
	*			 1 invalid.
	*/
	int TMeshToperator::flip22(HalfEdgeHandle &he_, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			/* check the halfedge handle is valid */
			if (he_.idx() < 0 || he_.idx() >= (int)tmesh->size_halfedge())
				throw std::range_error("FLIP2-2: The index of HalfFace handle is out of range!");

			// get two original tetras
			TetraHandle tetras[2];
			tetras[0] = tmesh->handle_to_entity(he_).hedron_handle();
			tetras[1] = tmesh->handle_to_entity(tmesh->radial_half_edge_handle(tmesh->mate_half_edge_handle(he_))).hedron_handle();

			// get the five related points
			PointHandle top, p1, p2, p[2];
			p1 = tmesh->point_handle(tmesh->from_vertex_handle(he_));
			p2 = tmesh->point_handle(tmesh->to_vertex_handle(he_));
			top = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(he_))));
			p[0] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(he_)));
			p[1] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(tmesh->mate_half_edge_handle(he_))))));

			// do flip
			// delete the two old tetrahedrons
			if (tmesh->erase_tetrahedron(tetras[0])==TetraHandle(-1))
				throw std::runtime_error("Tetrahedron erase fail!");
			if (tmesh->erase_tetrahedron(tetras[1])==TetraHandle(-1))
				throw std::runtime_error("Tetrahedron erase fail!");

			// add the two new tetrahedrons, the order of the points can be changed
			TetraHandle newTetra[2];
			newTetra[0] = tmesh->add_tetrahedron(tmesh->point(top), tmesh->point(p1), tmesh->point(p[1]), tmesh->point(p[0]));
			if (newTetra[0] == TetraHandle(-1))
				throw std::runtime_error("FLIP2-2: Tetrahedron add fail!");
			newTetra[1] = tmesh->add_tetrahedron(tmesh->point(top), tmesh->point(p2), tmesh->point(p[0]), tmesh->point(p[1]));
			if (newTetra[1] == TetraHandle(-1))
				throw std::runtime_error("FLIP2-2: Tetrahedron add fail!");

			// save recover data
			rec_flip22_top = top;
			rec_flip22_p[0] = p[0];
			rec_flip22_p[1] = p[1];
			rec_flip22_p1 = p1;
			rec_flip22_p2 = p2;
			rec_flip22_old_tetra[0] = tetras[0];
			rec_flip22_old_tetra[1] = tetras[1];
			rec_flip22_tetra[0] = newTetra[0];
			rec_flip22_tetra[1] = newTetra[1];

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra[0]);
				tetraVec_->push_back(newTetra[1]);
			}

		}
		catch (std::range_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;

	}

	/** topology operation : flip1-4.
	*	\param h_ is the index of the original tetra
	*   \param insertPoint is the coordinate of insert point
	*   \param tetraVec_(optional) return the new tetrahedrons
	*   \return  0 success.
	*			 1 invalid.
	*/
	int TMeshToperator::flip14(TetraHandle &h_, Point &insertPoint, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			/* check the halfedge handle is valid */
			if (h_.idx() < 0 || h_.idx() >= (int)tmesh->size_tetrahedron())
				throw std::range_error("FLIP1-4: The index of tetrahedron handle is out of range!");

			/* get points of the tetra*/
			PointHandle p[4];
			TetraMesh::HedronVertexIter hv_it;
			hv_it = tmesh->hedron_vertex_iter(h_);
			p[0] = tmesh->point_handle(hv_it.handle());++hv_it;
			p[1] = tmesh->point_handle(hv_it.handle());++hv_it;
			p[2] = tmesh->point_handle(hv_it.handle());++hv_it;
			p[3] = tmesh->point_handle(hv_it.handle());


			// do flip
			// delete the original tetrahedrons
			if (tmesh->erase_tetrahedron(h_)==TetraHandle(-1))
				throw std::runtime_error("Tetrahedron erase fail!");

			// add four new tetrahedrons, the order of the points can be changed
			TetraHandle newTetra[4];
			newTetra[0] = tmesh->add_tetrahedron(insertPoint, tmesh->point(p[0]), tmesh->point(p[2]), tmesh->point(p[1]));
			if (newTetra[0] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");
			newTetra[1] = tmesh->add_tetrahedron(insertPoint, tmesh->point(p[0]), tmesh->point(p[1]), tmesh->point(p[3]));
			if (newTetra[1] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");
			newTetra[2] = tmesh->add_tetrahedron(insertPoint, tmesh->point(p[0]), tmesh->point(p[3]), tmesh->point(p[2]));
			if (newTetra[2] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");
			newTetra[3] = tmesh->add_tetrahedron(insertPoint, tmesh->point(p[1]), tmesh->point(p[2]), tmesh->point(p[3]));
			if (newTetra[3] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");

			// save recover data
			rec_flip14_p[0] = p[0];
			rec_flip14_p[1] = p[1];
			rec_flip14_p[2] = p[2];
			rec_flip14_p[3] = p[3];
			rec_flip14_old_tetra = h_;
			rec_flip14_tetra[0] = newTetra[0];
			rec_flip14_tetra[1] = newTetra[1];
			rec_flip14_tetra[2] = newTetra[2];
			rec_flip14_tetra[3] = newTetra[3];

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra[0]);
				tetraVec_->push_back(newTetra[1]);
				tetraVec_->push_back(newTetra[2]);
				tetraVec_->push_back(newTetra[3]);
			}

		}
		catch (std::range_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;

	}

	
	/** topology operation : flip4-1.
	*	\param ph_ is the index of middle point
	*   \param tetraVec_(optional) return the new tetrahedrons
	*   \return  0 success.
	*			 1 invalid.
	*/
	int TMeshToperator::flip41(PointHandle &ph_, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			/* check the halfedge handle is valid */
			if (ph_.idx() < 0 || ph_.idx() >= (int)tmesh->size_point())
				throw std::range_error("FLIP4-1: The index of point handle is out of range!");

			/* get points of the tetra*/
			Point midpoint = tmesh->point(ph_);
			PointHandle p[4];                 /* the other four points */
			std::vector<TetraHandle> tetras;  /* tetras incident to the point*/
			std::set<PointHandle> pset;
			std::set<PointHandle>::iterator pset_it;
			TetraMesh::HedronVertexIter hv_it;
			tmesh->point_star(ph_, tetras);   /* get incident tetras */
			/* collecting point handles except ph_ */
			for (unsigned int i = 0; i < tetras.size(); ++i)
			{
				hv_it = tmesh->hedron_vertex_iter(tetras[i]);
				for (; hv_it; ++hv_it)
				{
					if (tmesh->point_handle(hv_it.handle()) != ph_)
						pset.insert(tmesh->point_handle(hv_it.handle()));
				}
			}
			pset_it = pset.begin();
			p[0] = *pset_it; ++pset_it;
			p[1] = *pset_it; ++pset_it;
			p[2] = *pset_it; ++pset_it;
			p[3] = *pset_it;

			// do flip
			// delete the original tetrahedrons
			for (unsigned int i = 0; i < tetras.size(); ++i)
			{
				if (tmesh->erase_tetrahedron(tetras[i])==TetraHandle(-1))
					throw std::runtime_error("Tetrahedron erase fail!");
			}

			// add a new tetrahedrons, the order of the points can be changed
			TetraHandle newTetra;
			newTetra = tmesh->add_tetrahedron(tmesh->point(p[0]), tmesh->point(p[1]), tmesh->point(p[2]), tmesh->point(p[3]));
			if (newTetra == TetraHandle(-1))
				throw std::runtime_error("FLIP4-1: Tetrahedron add fail!");

			// save recover data
			rec_flip41_midpoint = midpoint;
			rec_flip41_tetra = newTetra;
			rec_flip41_old_tetra[0] = tetras[0];
			rec_flip41_old_tetra[1] = tetras[1];
			rec_flip41_old_tetra[2] = tetras[2];
			rec_flip41_old_tetra[3] = tetras[3];

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra);
			}

		}
		catch (std::range_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	/** topology operation : flip1-2.
	*   \param he_ is the index of split edge
	*	\param ph_ is the index of middle point
	*   \param tetraVec_(optional) return the new tetrahedrons
	*   \return  0 success.
	*			 1 invalid.
	*/
	int TMeshToperator::flip12(HalfEdgeHandle &he_, PointHandle &ph_, std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			/* check the halfedge handle is valid */
			if (he_.idx() < 0 || he_.idx() >= (int)tmesh->size_halfedge())
				throw std::range_error("FLIP1-2: The index of HalfEdge handle is out of range!");
			if (ph_.idx() < 0 || ph_.idx() >= (int)tmesh->size_point())
				throw std::range_error("FLIP1-2: The index of point handle is out of range!");

			/* get points of the tetra*/
			PointHandle p[4];
			p[0] = tmesh->point_handle(tmesh->from_vertex_handle(he_));
			p[1] = tmesh->point_handle(tmesh->to_vertex_handle(he_));
			p[2] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(he_)));
			p[3] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(he_))));

			// do flip
			// delete the original tetrahedrons
			TetraHandle tetra = tmesh->handle_to_entity(he_).hedron_handle();
			if (tmesh->erase_tetrahedron(tetra) == TetraHandle(-1))
				throw std::runtime_error("Tetrahedron erase fail!");

			// add four new tetrahedrons
			TetraHandle newTetra[2];
			newTetra[0] = tmesh->add_tetrahedron(tmesh->point(ph_), tmesh->point(p[1]), tmesh->point(p[2]), tmesh->point(p[3]));
			if (newTetra[0] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");
			newTetra[1] = tmesh->add_tetrahedron(tmesh->point(ph_), tmesh->point(p[0]), tmesh->point(p[3]), tmesh->point(p[2]));
			if (newTetra[1] == TetraHandle(-1))
				throw std::runtime_error("FLIP3-2: Tetrahedron add fail!");

			// save recover data
			rec_flip12_p[0] = p[0];
			rec_flip12_p[1] = p[1];
			rec_flip12_p[2] = p[2];
			rec_flip12_p[3] = p[3];
			rec_flip12_old_tetra = tetra;
			rec_flip12_tetra[0] = newTetra[0];
			rec_flip12_tetra[1] = newTetra[1];

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(newTetra[0]);
				tetraVec_->push_back(newTetra[1]);
			}
		}
		catch (std::range_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	//********************************************END OF FLIP OPERATION*********************************************//

	//******************************************EDGE CONTRACTION OPERATION*******************************************//

	int TMeshToperator::edge_contract(HalfEdgeHandle &he_, PointHandle *newpoint)
	{
		//invalid handle
		if ( he_.idx() < 0 || he_.idx() >= (int)tmesh->size_halfedge())
			return 0;

		VertexHandle vf, vt;
		Point midpoint, pf_, pt_, pf, pt;
		PointHandle pht, phf; 
		std::vector<TetraHandle> edgeStarHH;
		std::vector<TetraHandle> vertexStar;
		std::vector<TetraHandle> tvstar;
		TetraMesh::HedronHalfEdgeIter hhe_it;
		Point tp[4];
		PointHandle ph_new;
		// get midpoint and endpoints of the edge
		phf = tmesh->point_handle(tmesh->from_vertex_handle(he_));
		pht = tmesh->point_handle(tmesh->to_vertex_handle(he_));
		pf_ = tmesh->point(tmesh->from_vertex_handle(he_));
		pt_ = tmesh->point(tmesh->to_vertex_handle(he_));
		midpoint = (pf_ + pt_)/ 2.0;
		// get tetras which containing the edge
		tmesh->edge_star(he_, edgeStarHH);
		// get tetras which containing the to_point
		tmesh->vertex_star(tmesh->to_vertex_handle(he_), tvstar);

		/* fetch tetras only incident to from_vertex */
		for (unsigned int i = 0; i < tvstar.size(); i++)
		{
			if (std::find(edgeStarHH.begin(), edgeStarHH.end(), tvstar[i]) == edgeStarHH.end())
				vertexStar.push_back(tvstar[i]);
		}

		// recovering data saving
		//rec_ec_valid = true;
		rec_ec_point = pf_;
		rec_ec_erase_point = pt_;
		rec_ec_point_handle = phf;
		rec_ec_to_point_vertex_container.clear();
		rec_ec_tetravertex.clear();
		rec_ec_old_tetra.clear();
		std::copy(edgeStarHH.begin(), edgeStarHH.end(), back_inserter(rec_ec_old_tetra));
		/* For each tetra saving two points, which are neither the from_point nor the to_point, 
		   as the recover data of the tetra. Because the from_point and to_point have been saved. */
		for (unsigned int i = 0; i < edgeStarHH.size(); i++)
		{
			TetraMesh::HedronVertexIter hv_it;
			PointHandle tmph;
			for (hv_it = tmesh->hedron_vertex_iter(edgeStarHH[i]); hv_it; ++hv_it)
			{
				tmph = tmesh->point_handle(hv_it.handle());
				if (tmph != phf && tmph != pht)
					rec_ec_tetravertex.push_back(tmph);
			}
		}

		//std::vector<VertexHandle>::iterator vec_it;
		//for (unsigned int i = 0; i < edgeStarHH.size(); i++)
		//{
		//	// find the contracting halfedge in the tetra & get endpoints
		//	for (hhe_it = tmesh->hedron_half_edge_iter(edgeStarHH[i]); hhe_it; ++hhe_it)
		//	{
		//		vf = tmesh->from_vertex_handle(hhe_it.handle());
		//		vt = tmesh->to_vertex_handle(hhe_it.handle());
		//		pf = tmesh->point(vf);
		//		pt = tmesh->point(vt);
		//		if (pf == pf_ && pt == pt_)
		//			break;
		//	}
		//	if (pf != pf_ || pt != pt_)
		//		return 0;

		//	// get halfFaces opposite to the endpoints of the contracting edge
		//	HalfFaceTH hff = tmesh->handle_to_entity(HalfFaceHandle(tmesh->from_vertex_handle(hhe_it.handle()).idx()));
		//	HalfFaceTH hft = tmesh->handle_to_entity(HalfFaceHandle(tmesh->to_vertex_handle(hhe_it.handle()).idx()));
		//	// get opposite halfFaces of hff and hft
		//	HalfFaceHandle ophff = hff.opposite_face_handle();
		//	HalfFaceHandle ophft = hft.opposite_face_handle();

		//	// get halfface pointhandle string which is used to build opposite halfface
		//	char buf1[128], buf2[128];
		//	PointHandle pointHandle[3];
		//	HalfEdgeTH he_e = tmesh->handle_to_entity(hff.first_half_edge_handle());
		//	pointHandle[0] = tmesh->start_point_handle(he_e.handle());
		//	pointHandle[1] = tmesh->start_point_handle(he_e.next_half_edge_handle());
		//	pointHandle[2] = tmesh->start_point_handle(he_e.prev_half_edge_handle());
		//	std::sort(pointHandle, pointHandle+3);
		//	sprintf_s(buf1,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);

		//	he_e = tmesh->handle_to_entity(hft.first_half_edge_handle());
		//	pointHandle[0] = tmesh->start_point_handle(he_e.handle());
		//	pointHandle[1] = tmesh->start_point_handle(he_e.next_half_edge_handle());
		//	pointHandle[2] = tmesh->start_point_handle(he_e.prev_half_edge_handle());
		//	std::sort(pointHandle, pointHandle+3);
		//	sprintf_s(buf2,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);

		//	// erase tetra
		//	tmesh->erase_tetrahedron(edgeStarHH[i]);

		//	// modify topology(reset opposite halfFace)
		//	if (ophff.idx() != -1)
		//		tmesh->handle_to_entity(ophff).set_opposite_face_handle(ophft);
		//	if (ophft.idx() != -1)
		//		tmesh->handle_to_entity(ophft).set_opposite_face_handle(ophff);

		//	// update opposite halfFace flag set _opphf[]
		//	tmesh->erase_opposite_face_connection(buf1);
		//	if (ophff.idx() == -1 && ophft.idx() == -1)
		//	{
		//		tmesh->erase_opposite_face_connection(buf2);
		//	}
		//	else if (ophff.idx() != -1 && ophft.idx() == -1)
		//	{
		//		tmesh->reset_oppsite_half_face(buf2, ophff.idx());
		//	}
		//	else if (ophff.idx() == -1 && ophft.idx() != -1)
		//	{
		//		tmesh->reset_oppsite_half_face(buf2, ophft.idx());
		//	}
		//	else
		//		tmesh->reset_oppsite_half_face(buf2, -1);
		//}


		// adjusting topology relation in the TetraMesh
		for (unsigned int i = 0; i <edgeStarHH.size(); i++)
		{
			tmesh->erase_tetrahedron(edgeStarHH[i]);
		}

		//----------------------edge contracting-------------------------//
		// move from_point to the midpoint
		ph_new = phf;
		tmesh->set_point(ph_new, midpoint);

		/* erase tetras incident to to_point
		   than add new tetras incident to from_point */
		TetraMesh::HedronVertexIter hv_it;
		TetraHandle newtet;
		int pcnt;
		/* fetch the first point of these new tetras */
		tp[0] = midpoint;
		for (unsigned int i = 0; i < vertexStar.size(); i++)
		{
			if (!tmesh->is_valid(vertexStar[i]))
				continue;

			pcnt = 1;
			/* fetch  */
			for (hv_it = tmesh->hedron_vertex_iter(vertexStar[i]); hv_it; ++hv_it)
			{	
				if (tmesh->point(hv_it.handle()) == pt_)
					continue;
				tp[pcnt++] = tmesh->point(hv_it.handle());
			}

			/* erase old one */
			tmesh->erase_tetrahedron(vertexStar[i]);

			/* add new one */
			newtet = tmesh->add_tetrahedron(tp[0], tp[1], tp[2], tp[3], &vertexStar[i]);

			/* save the left tetras which incident to the to_point */
			rec_ec_to_point_vertex_container.push_back(newtet);
		}

		// delete to_point
		tmesh->erase_point(pht);

		//------------------- end edge contracting------------------------//
		if (newpoint)
			*newpoint = ph_new;
		return 1;
	}

	//***************************************END OF EDGE CONTRACTION OPERATION***************************************//


	/************************************************************************/
	/* Tetrahedral Mesh Topological Operation                               */
	/* Recover Operation                                                    */
	/************************************************************************/

	/** operation flip1-2 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip12(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::runtime_error("No Mesh!");

			if (!rec_flip12_valid)
				throw std::runtime_error("RECOVER_FLIP12: Data not valid!");
			rec_flip12_valid = false;

			// delete new three tetrahedrons
			for (int i = 0; i < 2; ++i)
			{
				if (rec_flip12_tetra[i] == TetraHandle(-1))
					continue;
				if  (tmesh->erase_tetrahedron(rec_flip12_tetra[i]) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP12: Tetrahedron erase fail!");
			}

			// add the original tetrahedron
			TetraHandle tetra;
			tetra = tmesh->add_tetrahedron(tmesh->point(rec_flip12_p[0]), tmesh->point(rec_flip12_p[1]), 
				                           tmesh->point(rec_flip12_p[2]), tmesh->point(rec_flip12_p[3]),
										   &rec_flip12_old_tetra);
			if (tetra == TetraHandle(-1))
				throw std::runtime_error("RECOVER FLIP12: Tetrahedron add error!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetra);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
		}
		return 0;
	}

	/** operation flip1-3 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip13(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::runtime_error("No Mesh!");

			if (!rec_flip13_valid)
				throw std::runtime_error("RECOVER_FLIP13: Data not valid!");
			rec_flip13_valid = false;

			// delete new three tetrahedrons
			for (int i = 0; i < 3; ++i)
			{
				if (rec_flip13_tetra[i] == TetraHandle(-1))
					continue;
				if  (tmesh->erase_tetrahedron(rec_flip13_tetra[i]) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP13: Tetrahedron erase fail!");
			}

			// add the original tetrahedron
			TetraHandle tetra;
			tetra = tmesh->add_tetrahedron(tmesh->point(rec_flip13_p[0]), tmesh->point(rec_flip13_p[1]), 
				                           tmesh->point(rec_flip13_p[2]), tmesh->point(rec_flip13_p[3]),
										   &rec_flip13_old_tetra);
			if (tetra == TetraHandle(-1))
				throw std::runtime_error("RECOVER FLIP13: Tetrahedron add error!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetra);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
		}
		return 0;
	}

	/** operation flip2-3 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip23(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (!rec_flip23_valid)
				throw std::runtime_error("RECOVER_FLIP23: Data not valid!");
			rec_flip23_valid = false;

			// delete new three tetrahedrons
			for (int i = 0; i < 3; ++i)
			{
				if (rec_flip23_tetra[i] == TetraHandle(-1))
					continue;
				if  (tmesh->erase_tetrahedron(rec_flip23_tetra[i]) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP23: Tetrahedron erase fail!");
			}

			// add the original two tetrahedrons
			TetraHandle tetras[2];
			tetras[0] = tmesh->add_tetrahedron(tmesh->point(rec_flip23_top), tmesh->point(rec_flip23_p[0]), 
				                               tmesh->point(rec_flip23_p[1]), tmesh->point(rec_flip23_p[2]),
											   &rec_flip23_old_tetra[0]);
			if (tetras[0] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP23: Tetrahedron add fail!");

			tetras[1] = tmesh->add_tetrahedron(tmesh->point(rec_flip23_bot), tmesh->point(rec_flip23_p[2]), 
				                               tmesh->point(rec_flip23_p[1]), tmesh->point(rec_flip23_p[0]),
				                               &rec_flip23_old_tetra[1]);
			if (tetras[1] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP23: Tetrahedron add fail!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetras[0]);
				tetraVec_->push_back(tetras[1]);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}


	/** operation flip2-2 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip22(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (!rec_flip22_valid)
				throw std::runtime_error("RECOVER_FLIP22: Data not valid!");
			rec_flip22_valid = false;

			// delete new two tetrahedrons
			for (int i = 0; i < 2; ++i)
			{
				if (rec_flip22_tetra[i] == TetraHandle(-1))
					continue;
				if  (tmesh->erase_tetrahedron(rec_flip22_tetra[i]) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP22: Tetrahedron erase fail!");
			}

			// add the original two tetrahedrons
			TetraHandle tetras[2];
			tetras[0] = tmesh->add_tetrahedron(tmesh->point(rec_flip22_top), tmesh->point(rec_flip22_p1), 
				                               tmesh->point(rec_flip22_p[1]), tmesh->point(rec_flip22_p2),
											   &rec_flip22_old_tetra[0]);
			if (tetras[0] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP22: Tetrahedron add fail!");

			tetras[1] = tmesh->add_tetrahedron(tmesh->point(rec_flip22_top), tmesh->point(rec_flip22_p1), 
				                               tmesh->point(rec_flip22_p2), tmesh->point(rec_flip22_p[0]),
											   &rec_flip22_old_tetra[1]);
			if (tetras[1] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP22: Tetrahedron add fail!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetras[0]);
				tetraVec_->push_back(tetras[1]);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	/** operation flip3-2 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip32(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (!rec_flip32_valid)
				throw std::runtime_error("RECOVER_FLIP32: Data not valid!");
			rec_flip32_valid = false;

			// delete new two tetrahedrons
			for (int i = 0; i < 2; ++i)
			{
				if (rec_flip32_tetra[i] == TetraHandle(-1))
					continue;
				if  (tmesh->erase_tetrahedron(rec_flip32_tetra[i]) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP32: Tetrahedron erase fail!");
			}

			// add the original three tetrahedrons
			TetraHandle tetras[3];
			tetras[0] = tmesh->add_tetrahedron(tmesh->point(rec_flip32_top), tmesh->point(rec_flip32_p[0]), 
				                               tmesh->point(rec_flip32_p[1]), tmesh->point(rec_flip32_bot),
											   &rec_flip32_old_tetra[0]);
			if (tetras[0] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP32: Tetrahedron add fail!");

			tetras[1] = tmesh->add_tetrahedron(tmesh->point(rec_flip32_top), tmesh->point(rec_flip32_p[0]), 
				                               tmesh->point(rec_flip32_bot), tmesh->point(rec_flip32_p[2]),
											   &rec_flip32_old_tetra[1]);
			if (tetras[1] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP32: Tetrahedron add fail!");

			tetras[2] = tmesh->add_tetrahedron(tmesh->point(rec_flip32_top), tmesh->point(rec_flip32_bot), 
				                               tmesh->point(rec_flip32_p[1]), tmesh->point(rec_flip32_p[2]),
											   &rec_flip32_old_tetra[2]);
			if (tetras[2] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP32: Tetrahedron add fail!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetras[0]);
				tetraVec_->push_back(tetras[1]);
				tetraVec_->push_back(tetras[2]);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}


	/** operation flip3-1 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip31(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (!rec_flip31_valid)
				throw std::runtime_error("RECOVER_FLIP31: Data not valid!");
			rec_flip31_valid = false;

			// delete the new tetrahedron
			if (rec_flip31_tetra != TetraHandle(-1))
				if  (tmesh->erase_tetrahedron(rec_flip31_tetra) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP3-1: Tetrahedron erase fail!");

			// add the original three tetrahedrons
			TetraHandle tetras[3];
			tetras[0] = tmesh->add_tetrahedron(tmesh->point(rec_flip31_top), tmesh->point(rec_flip31_face_p[0]), 
				                               tmesh->point(rec_flip31_face_p[1]), rec_flip31_midpoint,
											   &rec_flip31_old_tetra[0]);
			if (tetras[0] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP31: Tetrahedron add fail!");

			tetras[1] = tmesh->add_tetrahedron(tmesh->point(rec_flip31_top), tmesh->point(rec_flip31_face_p[2]), 
				                               tmesh->point(rec_flip31_face_p[0]), rec_flip31_midpoint,
											   &rec_flip31_old_tetra[1]);
			if (tetras[1] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP31: Tetrahedron add fail!");

			tetras[2] = tmesh->add_tetrahedron(tmesh->point(rec_flip31_top), tmesh->point(rec_flip31_face_p[1]), 
				                               tmesh->point(rec_flip31_face_p[2]), rec_flip31_midpoint,
											   &rec_flip31_old_tetra[2]);
			if (tetras[2] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP31: Tetrahedron add fail!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetras[0]);
				tetraVec_->push_back(tetras[1]);
				tetraVec_->push_back(tetras[2]);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	/** operation flip1-4 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip14(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::runtime_error("No Mesh!");

			if (!rec_flip14_valid)
				throw std::runtime_error("RECOVER_FLIP14: Data not valid!");
			rec_flip14_valid = false;

			// delete new four tetrahedrons
			for (int i = 0; i < 4; ++i)
			{
				if (rec_flip14_tetra[i] == TetraHandle(-1))
					continue;
				if  (tmesh->erase_tetrahedron(rec_flip14_tetra[i]) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP14: Tetrahedron erase fail!");
			}

			// add the original tetrahedron
			TetraHandle tetra;
			tetra = tmesh->add_tetrahedron(tmesh->point(rec_flip14_p[0]), tmesh->point(rec_flip14_p[1]), 
			                               tmesh->point(rec_flip14_p[2]), tmesh->point(rec_flip14_p[3]),
										   &rec_flip14_old_tetra);
			if (tetra == TetraHandle(-1))
				throw std::runtime_error("RECOVER FLIP14: Tetrahedron add error!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetra);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
		}
		return 0;
	}

	/** operation flip4-1 recover
	*   param  : tetraVec_(optional) return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_flip41(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (!rec_flip41_valid)
				throw std::runtime_error("RECOVER_FLIP41: Data not valid!");
			rec_flip41_valid = false;

			/* get points of the new tetra */
			Point p[4];
			TetraMesh::HedronVertexIter hv_it;
			hv_it = tmesh->hedron_vertex_iter(rec_flip41_tetra);
			p[0] = tmesh->point(hv_it.handle()); ++ hv_it;
			p[1] = tmesh->point(hv_it.handle()); ++ hv_it;
			p[2] = tmesh->point(hv_it.handle()); ++ hv_it;
			p[3] = tmesh->point(hv_it.handle());

			// delete the new tetrahedron
			if (rec_flip41_tetra != TetraHandle(-1))
				if  (tmesh->erase_tetrahedron(rec_flip41_tetra) == TetraHandle(-1))
					throw std::runtime_error("RECOVER_FLIP4-1: Tetrahedron erase fail!");

			// add the original three tetrahedrons
			TetraHandle tetras[4];
			tetras[0] = tmesh->add_tetrahedron(rec_flip41_midpoint, p[0], p[2], p[1], &rec_flip41_old_tetra[0]);
			if (tetras[0] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP41: Tetrahedron add fail!");

			tetras[1] = tmesh->add_tetrahedron(rec_flip41_midpoint, p[0], p[1], p[3], &rec_flip41_old_tetra[1]);
			if (tetras[1] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP41: Tetrahedron add fail!");

			tetras[2] = tmesh->add_tetrahedron(rec_flip41_midpoint, p[0], p[3], p[2], &rec_flip41_old_tetra[2]);
			if (tetras[2] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP41: Tetrahedron add fail!");

			tetras[3] = tmesh->add_tetrahedron(rec_flip41_midpoint, p[1], p[2], p[3], &rec_flip41_old_tetra[3]);
			if (tetras[4] == TetraHandle(-1))
				throw std::runtime_error("RECOVER_FLIP41: Tetrahedron add fail!");

			if (tetraVec_ != NULL)
			{
				tetraVec_->push_back(tetras[0]);
				tetraVec_->push_back(tetras[1]);
				tetraVec_->push_back(tetras[2]);
				tetraVec_->push_back(tetras[3]);
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	/** operation recover
	*   \param: tetraVec_(optional): return the original tetrahedrons
	*   return : 1 succeed  
	*            0 failed
	*/
	int TMeshToperator::recover_edge_contract(std::vector<TetraHandle> *tetraVec_)
	{
		try
		{
			if (!tmesh)
				throw std::range_error("No Mesh!");

			if (!rec_ec_valid)
				throw std::runtime_error("RECOVER_EDGE_CONTRACT: Data not valid!");

			rec_ec_valid = false;

			PointHandle newph;
			std::vector<TetraHandle> pointStar;
			TetraMesh::HedronVertexIter hv_it;
			Point tp[4];
			TetraHandle newtet;

			/* get tetras incident to the reserved point */
			tmesh->point_star(rec_ec_point_handle, pointStar);   

			/* recover points */
			/* reset the reserved point */
			tmesh->set_point(rec_ec_point_handle, rec_ec_point);
			/* add the erased point back*/
			newph = tmesh->add_point(rec_ec_erase_point);

			// recover point vertex map of erased point
			if (rec_ec_to_point_vertex_container.size())
			{
				/* Iterate all the tetras used to incident to the to_point but now incident to the from_point,
				   than erase them and add new tetras incident to to_point
				*/
				int pcnt;
				tp[0] = rec_ec_erase_point;
				for (unsigned int i = 0; i < rec_ec_to_point_vertex_container.size(); i++)
				{
					pcnt = 1;
					/* iterate all the vertex to find the one that connects to the from_point */
					for (hv_it = tmesh->hedron_vertex_iter(rec_ec_to_point_vertex_container[i]); hv_it; ++hv_it)
					{
						if (tmesh->point_handle(hv_it.handle()) == rec_ec_point_handle)
							continue;
						tp[pcnt++] = tmesh->point(hv_it.handle());
					}
					/* erase old tetra */
					tmesh->erase_tetrahedron(rec_ec_to_point_vertex_container[i]);

					/* add new one */
					newtet = tmesh->add_tetrahedron(tp[0], tp[1], tp[2], tp[3], &rec_ec_to_point_vertex_container[i]);
				}
			}

			// reset opposite half face to leave room for recovering tetras
			//std::set<std::string> facestr;
			//std::set<std::string>::iterator fs_it;
			//TetraMesh::HedronFaceIter hf_it;
			//char buf[128];
			//PointHandle pointHandle[3];

			/* get points information of faces whose opposite half face should be changed */
			//for (unsigned int i = 0; i < rec_ec_tetravertex.size(); i = i+2)
			//{
			//	pointHandle[0] = rec_ec_tetravertex[i];
			//	pointHandle[1] = rec_ec_tetravertex[i+1];
			//	pointHandle[2] = newph;
			//	std::sort(pointHandle, pointHandle+3);
			//	sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
			//	facestr.insert(buf);

			//	pointHandle[0] = rec_ec_tetravertex[i];
			//	pointHandle[1] = rec_ec_tetravertex[i+1];
			//	pointHandle[2] = rec_ec_point_handle;
			//	std::sort(pointHandle, pointHandle+3);
			//	sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);
			//	facestr.insert(buf);
			//}

			/* reset opposite halfface for these faces */
			//int cnt = 0;
			//for (unsigned int i = 0; i < pointStar.size(); ++i)
			//{
			//	// for each face in the tetra
			//	for (hf_it = tmesh->hedron_face_iter(pointStar[i]); hf_it; ++hf_it)
			//	{
			//		// get its points information
			//		HalfEdgeTH he_e = tmesh->handle_to_entity(tmesh->handle_to_entity(hf_it.handle()).first_half_edge_handle());
			//		pointHandle[0] = tmesh->start_point_handle(he_e.handle());
			//		pointHandle[1] = tmesh->start_point_handle(he_e.next_half_edge_handle());
			//		pointHandle[2] = tmesh->start_point_handle(he_e.prev_half_edge_handle());
			//		std::sort(pointHandle, pointHandle+3);
			//		sprintf_s(buf,"%d %d %d",pointHandle[0],pointHandle[1],pointHandle[2]);

			//		// check if the halfface needs to change its opposite halfface
			//		fs_it = find(facestr.begin(), facestr.end(), buf);
			//		if (fs_it != facestr.end())
			//		{
			//			tmesh->reset_oppsite_half_face(buf, hf_it.handle().idx());
			//			tmesh->handle_to_entity(hf_it.handle()).set_opposite_face_handle(HalfFaceHandle(-1));
			//		
			//			++ cnt;
			//			break;
			//		}
			//	}
			//}

			// recover tetrahedrons
			std::vector<TetraHandle> tetras;
			TetraHandle tetra;

			for (unsigned int i = 0; i < rec_ec_tetravertex.size(); i = i+2)
			{
				/* add tetras back to mesh 
				 * each tetra contains four point: the two endpoints of the edge(rec_ec_point and rec_ec_erase_point),
				 * and the other two points reserved in rec_ec_tetravertex */
				tetra = tmesh->add_tetrahedron(rec_ec_point, rec_ec_erase_point, 
					                           tmesh->point(rec_ec_tetravertex[i]), 
					                           tmesh->point(rec_ec_tetravertex[i+1]),
											   &rec_ec_old_tetra[i/2]);
				if (tetra == TetraHandle(-1))
					throw std::runtime_error("RECOVER_EC: Tetrahedron add fail!");
				tetras.push_back(tetra);
			}

			// set output tetras set 
			if (tetraVec_ != NULL)
			{
				copy(tetras.begin(), tetras.end(), std::back_inserter(*tetraVec_));
			}
		}
		catch (std::runtime_error& e)
		{
			std::cout << e.what() << std::endl;
			return 1;
		}
		return 0;
	}

	/************************************************************************/
	/* Topology pass                                                        */
	/************************************************************************/

	/* perform a pass of topological improvement.
	for now, this means trying to remove each edge
	of each tet in the stack that is passed in,
	and if no edge can be removed, trying to remove
	each face. */
	bool TMeshToperator::topologyPass(std::vector<TetraHandle> *tetstack, std::vector<TetraHandle> *outstack, double &minqualafter)
	{
		double outquality;                      /* the quality of the current tet */
		TetraHandle stacktet;                   /* stack tet */
		int removed, flipped;                   /* was edge / face removal successful */
		//int numremoved = 0;                     /* number of edges removed */
		//int numflipped = 0;                     /* number of 2-3 face flips */
		//int biggestring = 0;                    /* biggest ring seen */
		//int biggestface = 0;
		std::vector<TetraHandle> outtet;        /* stack tet : generated by single step operation */
		std::vector<TetraHandle> savestack;     /* save the input stack in case of failure */
		std::set<TetraHandle> incidenttet;   /* stack tet : influenced tetras in topology pass */
		double minqualbefore;
		int tetcnt;

		minqualafter = HUGEFLOAT;

		/* reset output stack */
		if (outstack != NULL)
		{
			outstack->clear();
		}

		/* compute worst input quality.*/
		if (tetstack != NULL)
		{
			/* save a copy of the input stack */
			copy(tetstack->begin(), tetstack->end(), back_inserter(savestack));
			minqualbefore = qualcalculator->minstackquality(tmesh, *tetstack, qualityMetric);
		}
		else
		{
			minqualbefore = qualcalculator->minmeshquality(tmesh, qualityMetric);
		}

		/* if we aren't allowed to do any topo stuff, don't even bother */
		if (edgeremoval == 0 &&
			singlefaceremoval == 0 &&
			multifaceremoval == 0)
		{
			minqualafter = minqualbefore;
			if (outstack != NULL)
			{
				outstack->clear();
				copy(savestack.begin(), savestack.end(), back_inserter(*outstack));
			}
			return false;
		}

		/* fetch the number of tetras */
		if (tetstack != NULL)
		{
			tetcnt = tetstack->size();
		}
		else
			tetcnt = tmesh->size_tetrahedron();

		/* topology improvement: go through each tet on the stack */
		for (int i = 0; i < tetcnt; ++i)
		{
			/* pull the top tet off the stack */
			if (tetstack != NULL)
				stacktet = (*tetstack)[i];
			else
				stacktet = TetraHandle(i);

			if (!tmesh->is_valid(stacktet))
				continue;

			if(i == 269)
				int a = 0;

			if (edgeremoval)
			{
				/* try edge removal first */
				outtet.clear();
				removed = tryedgeremove(stacktet, outtet, outquality);
				++stats->edgeremovalattempts;
			}
			else
			{
				removed = false;
			}

			if (removed)
				++stats->edgeremovals;

			/* if edge removal failed, try face removal */
			if (!removed)
			{
				/* now try multi-face removal */
				outtet.clear();
				flipped = tryremovefaces(stacktet, outtet, outquality);
				++stats->faceremovalattempts;
			}
			else
			{
				flipped = false;
			}

			if (flipped)
				++stats->faceremovals;

			///* if nothing was done to this tet, push it on the output stack */
			//if (!removed && !flipped)
			//{
			//	/* push this tet on the output stack */
			//	outtet.clear();
			//	outtet.push_back(stacktet);
			//}

			/* fetch the influenced tetras */
			//copy(outtet.begin(), outtet.end(), back_inserter(incidenttet));

			if (removed || flipped)
			{
				for (unsigned int idx = 0; idx < outtet.size(); idx++)
				{
					incidenttet.insert(outtet[idx]);
				}
			}
			else
				incidenttet.insert(stacktet);

			removed = false;
			flipped = false;

		}

		/* calculate the new worst quality of tetras */
		outtet.clear();
		copy(incidenttet.begin(), incidenttet.end(), back_inserter(outtet));
		if (tetstack != NULL)
			minqualafter = qualcalculator->minstackquality(tmesh, outtet, qualityMetric);
		else
			minqualafter = qualcalculator->minmeshquality(tmesh, qualityMetric);

		if (outstack != NULL)
			copy(outtet.begin(), outtet.end(), back_inserter(*outstack));

		/* if the new quality is better then set outstack and retrun true */
		if (minqualafter > minqualbefore)
			return true;

		return false;
	}

	/* attempt to remove the edge from vtx1 to vtx2. do flip22 or flip32.
	 * [PARAM] he_ : is the edge to be removed
	 * [PARAM] oldminqual : original worst quality
	 * [PARAM] newtets : new generated tetras
	 * [PARAM] outminqual : the worst quality after edge removal
	 * [PARAM] boundary : indicate if the edge is a boundary one
	 * [RETRUN] bool : TRUE edge removal success
	 *                 FALSE failed
	*/
	bool TMeshToperator::removeedge(HalfEdgeHandle he_, double &oldminqual, std::vector<TetraHandle> &newtets, 
			                        double &outminqual, bool &boundary)
	{
		double minqual = 1.0;  /*worst quality*/
		double tempqual;
		std::vector<TetraHandle> edgestar; /* tetras incident to he_*/
		tmesh->edge_star(he_, edgestar);
	
		/* do nothing if the condition does not satisfy our requirement */
		if (edgestar.size()!=2 && edgestar.size()!=3)
			return false;
	
		/* calculate the original worst quality */
		oldminqual = qualcalculator->minstackquality(tmesh, edgestar, qualityMetric);

		/* Before we actually do the flip, 
		   we should check the worst quality of new tetras 
		   is better than the original one's.*/
		/* do flip22 */
		if (edgestar.size() == 2 && tmesh->is_boundary(he_))
		{
			boundary = true;

			/*check if the faces beside the edge are coplanar */
			TetraMesh::HalfFaceHandle hf[2];
			hf[0] = tmesh->handle_to_entity(he_).half_face_handle();
			hf[1] = tmesh->mate_half_edge(he_).half_face_handle();
			if ((tmesh->normal(hf[0])-tmesh->normal(hf[1])).norm() > DEPSILON)
				return false;

			// get the five related point
			PointHandle top, p1, p2, p[2];  /* points of the tetras */
			p1 = tmesh->point_handle(tmesh->from_vertex_handle(he_));
			p2 = tmesh->point_handle(tmesh->to_vertex_handle(he_));
			top = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(he_))));
			p[0] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(he_)));
			p[1] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(tmesh->mate_half_edge_handle(he_))))));

			/* calculate the quality of one of the new tetra*/
			tempqual = qualcalculator->tetquality(tmesh->point(top), tmesh->point(p1), tmesh->point(p[1]), tmesh->point(p[0]), qualityMetric);
			if (tempqual < minqual)
				minqual = tempqual;

			/* calculate the quality of the other one of the new tetra*/
			tempqual = qualcalculator->tetquality(tmesh->point(top), tmesh->point(p2), tmesh->point(p[0]), tmesh->point(p[1]), qualityMetric);
			if (tempqual < minqual)
				minqual = tempqual;

			/* if the worst quality after flip is worse than the original one, than flip failed */
			if (minqual < oldminqual)
				return false;

			/* do flip22 */
			flip22(he_, &newtets);
			/* record this vertex insertion in the journal */
			jpointvec.clear();jthvec.clear();jphvec.clear();jvhvec.clear();
			getRecoverData(FLIP22, jpointvec, jthvec, jphvec, jvhvec);
			optjournal->pushJournal(FLIP22, jpointvec, jthvec, jphvec, jvhvec);

		}
		/* do flip23 */
		else
		{
			if(tmesh->is_boundary(he_))
				return false;

			boundary = false;

			// get the five related points
			PointHandle top, bot, p[3];  /* points of the tetras */
			top = tmesh->point_handle(tmesh->from_vertex_handle(he_));
			bot = tmesh->point_handle(tmesh->to_vertex_handle(he_));
			p[0] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(he_)));
			p[1] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(he_))));
			p[2] = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(he_)))));

			/* calculate the quality of the first tetra*/
			tempqual = qualcalculator->tetquality(tmesh->point(top), tmesh->point(p[0]), tmesh->point(p[1]), tmesh->point(p[2]), qualityMetric);
			if (tempqual < minqual)
				minqual = tempqual;

			/* calculate the quality of the second tetra*/
			tempqual = qualcalculator->tetquality(tmesh->point(bot), tmesh->point(p[2]), tmesh->point(p[1]), tmesh->point(p[0]), qualityMetric);
			if (tempqual < minqual)
				minqual = tempqual;

			/* if the worst quality after flip is worse than the original one, than flip failed */
			if (minqual < oldminqual)
				return false;

			/* do flip32 */
			flip32(he_, &newtets);
			/* record this vertex insertion in the journal */
			jpointvec.clear();jthvec.clear();jphvec.clear();jvhvec.clear();
			getRecoverData(FLIP32, jpointvec, jthvec, jphvec, jvhvec);
			optjournal->pushJournal(FLIP32, jpointvec, jthvec, jphvec, jvhvec);

		}

		outminqual = minqual;
		return true;
	}

	/* attempt to remove the face hf_. do flip23.
	[PARAM] hf_ : is the face to be removed
	[PARAM] oldminqual : original worst quality
	[PARAM] newtets : new generated tetras
	[PARAM] outminqual : the worst quality after edge removal
	[PARAM] boundary : indicate if the edge is a boundary one
	[RETRUN] bool : TRUE edge removal success
	FALSE failed
	*/
	bool TMeshToperator::removeface(HalfFaceHandle hf_, double oldminqual, std::vector<TetraHandle> &newtets,
		                            double &outminqual, bool &boundary)
	{
		double minqual = 1.0;  /*worst qualiyt*/
		double tempqual;
		bool isflip22 = false;
		HalfEdgeHandle he_;
		TetraMesh::FaceHalfedgeIter fe_iter;
		TetraMesh::Normal edgenorm;

		// get related five points
		PointHandle top, bot, p[3];
		fe_iter = tmesh->face_half_edge_iter(hf_);
		top = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(fe_iter.handle()))));
		bot = tmesh->point_handle(tmesh->to_vertex_handle(tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(tmesh->radial_half_edge_handle(fe_iter.handle())))));
		p[0] = tmesh->point_handle(tmesh->from_vertex_handle(fe_iter.handle()));++fe_iter;
		p[1] = tmesh->point_handle(tmesh->from_vertex_handle(fe_iter.handle()));++fe_iter;
		p[2] =tmesh-> point_handle(tmesh->from_vertex_handle(fe_iter.handle()));

		/* calculate the quality of first tetra*/
		tempqual = qualcalculator->tetquality(tmesh->point(top), tmesh->point(p[0]), tmesh->point(p[1]), tmesh->point(bot), qualityMetric);
		if (tempqual < minqual)
			minqual = tempqual;

		/* calculate the quality of second tetra*/
		tempqual = qualcalculator->tetquality(tmesh->point(top), tmesh->point(p[0]), tmesh->point(bot), tmesh->point(p[2]), qualityMetric);
		if (tempqual < minqual)
			minqual = tempqual;

		/* calculate the quality of third tetra*/
		tempqual = qualcalculator->tetquality(tmesh->point(top), tmesh->point(bot), tmesh->point(p[1]), tmesh->point(p[2]), qualityMetric);
		if (tempqual < minqual)
			minqual = tempqual;

		/* if the worst quality after flip is worse than the original one, than flip failed */
		if (minqual < oldminqual)
			return false;

		/* do flip22 */
		flip23(hf_, &newtets);

		/* record this vertex insertion in the journal */
		jpointvec.clear();jthvec.clear();jphvec.clear();jvhvec.clear();
		getRecoverData(FLIP23, jpointvec, jthvec, jphvec, jvhvec);
		optjournal->pushJournal(FLIP23, jpointvec, jthvec, jphvec, jvhvec);

		outminqual = minqual;
		return true;
	}


	/* try edge removal on all 6 edges of the specified tet. if removal
	succeeds, return 1. Otherwise, return 0. */
	bool TMeshToperator::tryedgeremove(TetraHandle tetrahandle, std::vector<TetraHandle> &newtetras, double &minafterqual)
	{
		/* first, check to make sure this tet still exists in the mesh */
		if (!tmesh->is_valid(tetrahandle))
			return false;

		bool edgeremoved = false;
		double oldminqual, afterminqual;
		bool boundary;
		std::vector<TetraHandle> newtets;
		std::vector<TetraHandle> edgestar;
		TetraMesh::HedronHalfEdgeIter he_iter;

		/* check each of this tet's six edges to see if it is eligible for 2-2 flip or 2-3 flip */
		for (he_iter = tmesh->hedron_half_edge_iter(tetrahandle); he_iter; ++he_iter)
		{
			/* don't try other edges if previous one has been removed */
			if (edgeremoved)
			{
				/* save data */
				newtetras = newtets;
				minafterqual = afterminqual;
				break;
			}
		
			/* do edge removal */
			newtets.clear();
			edgeremoved = removeedge(he_iter.handle(), oldminqual, newtets, afterminqual, boundary);
		}

		return edgeremoved;
	}


	/* try a 2-3 flip on all four faces of the specified tet
	returns 1 if a flip was performed, 0 otherwise */
	bool TMeshToperator::tryremovefaces(TetraHandle tetrahandle, std::vector<TetraHandle> &newtetras, double &minafterqual)
	{
		/* first, check to make sure this tet still exists in the mesh */
		if (!tmesh->is_valid(tetrahandle))
			return false;

		bool faceremoved = false;
		bool boundary;
		double oldminqual, afterminqual;
		TetraHandle opph;
		std::vector<TetraHandle> oldtets;
		std::vector<TetraHandle> newtets;
		TetraMesh::HedronFaceIter hf_iter;

		/* check each of this tet's four faces to see if it is eligible for 2-3 flip */
		for (hf_iter = tmesh->hedron_face_iter(tetrahandle); hf_iter; ++ hf_iter)
		{
			/* don't try other faces if previous one has been removed */
			if (faceremoved)
			{
				/* save data */
				newtetras = newtets;
				minafterqual = afterminqual;
				break;
			}

			/* get the opposite tetra of the face.
			   if the tetra is not exist then skip this face.*/
			if (!tmesh->has_opposite_half_face(hf_iter.handle()))
				continue;
			opph = tmesh->opposite_half_face(hf_iter.handle()).hedron_handle();
			if (!tmesh->is_valid(opph))
				continue;
			/* calculate the worst quality of old tetras */
			oldtets.clear();
			oldtets.push_back(tetrahandle);
			oldtets.push_back(opph);
			oldminqual = qualcalculator->minstackquality(tmesh, oldtets, qualityMetric);

			/* do face removal */
			newtets.clear();
			faceremoved = removeface(hf_iter.handle(), oldminqual, newtets, afterminqual, boundary);
		}

		return faceremoved;
	}

	/* go after the worst tets with contraction */
	void TMeshToperator::contractworst(double percentinsert, double bestmeans[], double outmeanqual[], double *outminqual, bool desperate)
	{
		double minqual;/* quality of the worst tet in the mesh */
		std::vector<TetraHandle> tetstack;          /* stack of tets  */
		double meshworstqual, origmeshworstqual;
		double fillthresh;

		qualcalculator->meshquality(tmesh, qualityMetric, &minqual);

		if (desperate == false)
		{    
			/* fill the stack of with the worst percent insert tets */
			fillstackpercent(tmesh, &tetstack, qualityMetric, 1.0);
		}
		else
		{
			if (minqual + QUALFROMDESPERATE < improvebehave->maxinsertquality[improvebehave->qualmeasure])
			{
				fillthresh = minqual + QUALFROMDESPERATE;
			}
			else if (minqual > improvebehave->maxinsertquality[improvebehave->qualmeasure])
			{
				fillthresh = minqual + QUALUPFROMDESPERATE;
			}
			else
			{
				fillthresh = improvebehave->maxinsertquality[improvebehave->qualmeasure];
			}
			fillstackqual(tmesh, &tetstack, qualityMetric, fillthresh);
		}

		origmeshworstqual = meshworstqual = minqual;

		/* perform insertion pass on stack of tets */
		contractpass(&tetstack, NULL, *outminqual, false);

		qualcalculator->meshquality(tmesh, qualityMetric, outminqual);
	}

	/* given a mesh and a percentage p,
	return the worst numtets * p tets in the mesh */
	void TMeshToperator::fillstackpercent(TetraMesh *mesh, std::vector<TetraHandle> *stack, int qualmeasure, double percent)
	{
		Point p[4];
		TetraMesh::HedronIter h_iter;
		int numouttets = 0;

		/* create a local stack for the global build */
		TetraStack alltetstack;
		ImproveTetra tetras;

		assert(percent > 0.0 && percent <= 1.0);

		/* fill it with every tet in the mesh */
		for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++ h_iter)
		{
			tetPoints(mesh, h_iter.handle(), p);
			tetras.handle = h_iter.handle();
			tetras.quality = qualcalculator->tetquality(p[0], p[1], p[2], p[3], qualityMetric);
			alltetstack.push_back(tetras);
		}

		/* sort it from worst to best */
		sort(alltetstack.begin(), alltetstack.end());

		/* figure out how many tets to take */
		numouttets = (int) (percent * ((double) alltetstack.size() + 1.0));

		if (numouttets > MAXINSERTTETS)
		{
			numouttets = MAXINSERTTETS;
		}

		/* make sure the output stack is empty */
		stack->clear();

		/* now, pop numouttets off the output stack */
		for (int i = 0; i < numouttets; i++)
		{
			stack->push_back(alltetstack[i].handle);
		}
	}

	void TMeshToperator::fillstackqual(TetraMesh *mesh, std::vector<TetraHandle> *stack, int qualmeasure,
		                               double threshold)
	{
		Point p[4];
		TetraMesh::HedronIter h_iter;
		int numouttets = 0;

		/* create a local stack for the global build */
		TetraStack alltetstack;
		ImproveTetra tetras;

		/* fill it with every tet in the mesh */
		for (h_iter = mesh->hedrons_begin(); h_iter != mesh->hedrons_end(); ++ h_iter)
		{
			tetPoints(mesh, h_iter.handle(), p);
			tetras.handle = h_iter.handle();
			tetras.quality = qualcalculator->tetquality(p[0], p[1], p[2], p[3], qualityMetric);
			alltetstack.push_back(tetras);
		}

		/* sort it from worst to best */
		sort(alltetstack.begin(), alltetstack.end());

		if (numouttets > MAXINSERTTETS)
		{
			numouttets = MAXINSERTTETS;
		}

		/* make sure the output stack is empty */
		stack->clear();

		/* now, pop numouttets off the output stack */
		for (int i = 0; alltetstack[i].quality < threshold; i++)
		{
			stack->push_back(alltetstack[i].handle);
		}
	}


	/* for each tet in the stack, try to contract its edges */
	bool TMeshToperator::contractpass(std::vector<TetraHandle> *tetstack, std::vector<TetraHandle> *outstack, 
		                              double &minqualafter, bool justfirstedge)
	{
		double outquality;           /* the quality of the current tet */
		TetraHandle stacktet;   /* point to stack tet */
		std::vector<TetraHandle> outtet;     /* tets influenced in one single step */
		//std::vector<TetraHandle> savestack;
		std::vector<TetraHandle> incidenttet;
		std::vector<TetraHandle>::iterator tet_iter;
		bool contracted;                /* was edge / face removal successful */
		//int origstacksize;             /* number of tets in the original stack */
		//int beforeid = optjournal->size();
		double minqualbefore, meanqualbefore[NUMMEANTHRESHOLDS];
		double minnow;
		TetraMesh::HedronVertexIter hv_it;
		bool validtet;

		minqualafter = HUGEFLOAT;

		//origstacksize = tetstacks + 1;

		/* compute worst input quality. do it global if no output stack */
		if (tetstack != NULL)
		{
			minqualbefore = qualcalculator->minstackquality(tmesh, *tetstack, qualityMetric);
			/* save a copy of the input stack */
			//copy(tetstack->begin(), tetstack->end(), back_inserter(savestack));
		}
		else
		{
			minqualbefore = qualcalculator->minmeshquality(tmesh, qualityMetric);
		}

		/* reset output stack */
		if (outstack != NULL)
		{
			outstack->clear();
		}

		/* fetch the number of tetras */
		int tetcnt;
		if (tetstack!=NULL)
			tetcnt = tetstack->size();
		else
			tetcnt = tmesh->size_tetrahedron();

		/* do edge contraction : go through each tet on the stack */
		for (int tetidx = 0; tetidx < tetcnt; ++tetidx)
		{
			/* pull the top tet off the stack */
			if (tetstack != NULL)
				stacktet = (*tetstack)[tetidx];
			else
				stacktet = TetraHandle(tetidx);

			if (!tmesh->is_valid(stacktet))
				continue;

			/* check if the tetra has boundary point,
			   if it has, then skip this tetra */
			validtet = true;
			for (hv_it = tmesh->hedron_vertex_iter(stacktet); hv_it; ++hv_it)
			{
				if (tmesh->is_boundary(hv_it.handle()))
				{
					validtet = false;
					break;
				}
			}
			if (!validtet)
				continue;

			/* try to contract this edge */
			outtet.clear();
			contracted = tryedgecontract(stacktet, &outtet, &outquality);

			if (contracted) 
				++stats->edgecontractions;

			/* if nothing was done to this tet, push it on the output stack */
			if (!contracted)
			{
				/* push this tet on the output stack */
				outtet.push_back(stacktet);
			}

			/* fetch the influenced tetras */
			copy(outtet.begin(), outtet.end(), back_inserter(incidenttet));

			contracted = false;
		}

		/* calculate quality after edge contraction */
		if (tetstack != NULL)
			minqualafter = qualcalculator->minstackquality(tmesh, incidenttet, qualityMetric);
		else
			minqualafter = qualcalculator->minmeshquality(tmesh, qualityMetric);

		/* fetch the output tetras */
		if (outstack != NULL)
		{
			copy(incidenttet.begin(), incidenttet.end(), back_inserter(*outstack));
		}

		/* if we didn't improve, undo this pass and indicate failure */
		if (minqualafter - minqualbefore < MINCONTRACTIMPROVEMENT)
			return false;

		return true;
	}


	/* try to contract all six edges of a tet */
	bool TMeshToperator::tryedgecontract(TetraHandle tetrahandle, std::vector<TetraHandle> *outstack, double *minqualafter)
	{
		Point p[4];
		bool contracted = false;
		double minqualbefore, meanqualbefore[NUMMEANTHRESHOLDS];
		double minqualnew;
		PointHandle newpoint;
		Point contractedge[6][2];
		//TetraMesh::HalfEdgeHandle contractedge[6];
		std::vector<TetraHandle> stacktet;
		std::vector<TetraHandle> pointstar;
		std::vector<TetraHandle> edgestar;
		TetraMesh::HalfEdgeHandle heh;
		TetraMesh::HedronHalfEdgeIter he_it;
		Point ep[2];
		Record record;
		Point pf, pt, pmid;
		std::vector<TetraHandle> tpstar;
		std::vector<TetraHandle> fpstar;
		std::vector<TetraHandle> newstar;
		TetraMesh::HedronVertexIter hv_it;
		TetraMesh oldmesh;
		oldmesh = *tmesh;

		//testing
		double minmeshqual;
		minmeshqual = qualcalculator->minmeshquality(tmesh, qualityMetric);

		/* first, check to make sure this tet still exists in the mesh */
		if (!tmesh->is_valid(tetrahandle))
		{
			return false;
		}

		/* get half edges need to be contracted */

		//contractedge[0] = TetraMesh::HalfEdgeHandle(tetrahandle.idx()*12);
		//contractedge[1] = tmesh->next_half_edge_handle(contractedge[0]);
		//contractedge[2] = tmesh->next_half_edge_handle(contractedge[1]);
		//contractedge[3] = tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(contractedge[0]));
		//contractedge[4] = tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(contractedge[1]));
		//contractedge[5] = tmesh->next_half_edge_handle(tmesh->mate_half_edge_handle(contractedge[2]));

		/* loop over all 6 edges of the tetrahedron */
		contractedge[0][0] = tmesh->point(VertexHandle(tetrahandle.idx()*4));
		contractedge[0][1] = tmesh->point(VertexHandle(tetrahandle.idx()*4+1));
		contractedge[1][0] = tmesh->point(VertexHandle(tetrahandle.idx()*4));
		contractedge[1][1] = tmesh->point(VertexHandle(tetrahandle.idx()*4+2));
		contractedge[2][0] = tmesh->point(VertexHandle(tetrahandle.idx()*4));
		contractedge[2][1] = tmesh->point(VertexHandle(tetrahandle.idx()*4+3));
		contractedge[3][0] = tmesh->point(VertexHandle(tetrahandle.idx()*4+1));
		contractedge[3][1] = tmesh->point(VertexHandle(tetrahandle.idx()*4+2));
		contractedge[4][0] = tmesh->point(VertexHandle(tetrahandle.idx()*4+1));
		contractedge[4][1] = tmesh->point(VertexHandle(tetrahandle.idx()*4+3));
		contractedge[5][0] = tmesh->point(VertexHandle(tetrahandle.idx()*4+2));
		contractedge[5][1] = tmesh->point(VertexHandle(tetrahandle.idx()*4+3));

		/* loop over all 6 edges of the tetrahedron */
		for (int idx = 0; idx < 6; idx++) 
		{
			contracted = true;
			minqualnew = 1.0;

			/* fetch current halfedge */
			heh = HalfEdgeHandle(-1);
			for (he_it = tmesh->hedron_half_edge_iter(tetrahandle); he_it; ++ he_it)
			{
				ep[0] = tmesh->point(tmesh->from_vertex_handle(he_it.handle()));
				ep[1] = tmesh->point(tmesh->to_vertex_handle(he_it.handle()));
				if ((ep[0] == contractedge[idx][0] && ep[1] == contractedge[idx][1]) || 
					(ep[1] == contractedge[idx][0] && ep[0] == contractedge[idx][1]))
				{
					heh = he_it.handle();
					break;
				}
			}
			//heh = contractedge[idx];

			if (heh.idx() == -1)
				return false;

			/* fetch tetras incident to edge */
			edgestar.clear();
			//star.clear();
			tmesh->edge_star(heh, edgestar);
			minqualbefore = qualcalculator->minstackquality(tmesh, edgestar, qualityMetric);

			//************************************************************************/
			//* check if edge contraction could improve quality                      */
			//************************************************************************/
			/* fetch from point */
			pf = tmesh->point(tmesh->from_vertex_handle(heh));
			/* fetch to point */
			pt = tmesh->point(tmesh->to_vertex_handle(heh));
			/* calculate the middle point */
			pmid = (pf + pt) / 2.0;

			/* fetch tetras incident to from_point */
			fpstar.clear();
			tmesh->vertex_star(tmesh->from_vertex_handle(heh), fpstar);
			/* fetch tetras incident to to_point */
			tpstar.clear();
			tmesh->vertex_star(tmesh->to_vertex_handle(heh), tpstar);

			/* get the original worst quality */
			//copy(fpstar.begin(), fpstar.end(), back_inserter(star));
			//copy(tpstar.begin(), tpstar.end(), back_inserter(star));

			/* move from_point and to_point to mid_point */
			tmesh->set_point(tmesh->point_handle(tmesh->from_vertex_handle(heh)), pmid);
			tmesh->set_point(tmesh->point_handle(tmesh->to_vertex_handle(heh)), pmid);

			/* fetch the tetras incident to from_point and to_point 
			   but do not incident to the contracting edge */
			newstar.clear();
			for (unsigned int i = 0; i < fpstar.size(); ++i)
			{
				if (std::find(edgestar.begin(), edgestar.end(), fpstar[i]) == edgestar.end())
				{
					newstar.push_back(fpstar[i]);
				}
			}
			for (unsigned int i = 0; i < tpstar.size(); ++i)
			{
				if (std::find(edgestar.begin(), edgestar.end(), tpstar[i]) == edgestar.end())
				{
					newstar.push_back(tpstar[i]);
				}
			}

			/* get new quality */
			minqualnew = qualcalculator->minstackquality(tmesh, fpstar, qualityMetric);

			minqualnew = qualcalculator->minstackquality(tmesh, tpstar, qualityMetric);

			minqualnew = qualcalculator->minstackquality(tmesh, newstar, qualityMetric);

			/* if the operator could not improve mesh quality, 
			   reset the two endpoints to their original positions
			   */
			if (minqualnew <= minqualbefore)
			{
				tmesh->set_point(tmesh->point_handle(tmesh->from_vertex_handle(heh)), pf);
				tmesh->set_point(tmesh->point_handle(tmesh->to_vertex_handle(heh)), pt);

				minmeshqual = qualcalculator->minmeshquality(tmesh, qualityMetric);

				minqualbefore = qualcalculator->minstackquality(tmesh, edgestar, qualityMetric);

				continue;
			}

			/* do edge contraction */
			tmesh->set_point(tmesh->point_handle(tmesh->from_vertex_handle(heh)), pf);
			tmesh->set_point(tmesh->point_handle(tmesh->to_vertex_handle(heh)), pt);

			minqualbefore = qualcalculator->minstackquality(tmesh, edgestar, qualityMetric);

			edge_contract(heh, &newpoint);

			tmesh->point_star(newpoint, edgestar);

			minqualbefore = qualcalculator->minstackquality(tmesh, edgestar, qualityMetric);

			minmeshqual = qualcalculator->minmeshquality(tmesh, qualityMetric);

			return true;



			///* if edge contraction succeeded */
			//if (edge_contract(heh, &newpoint))
			//{
			//	/* record this vertex insertion in the journal */
			//	jpointvec.clear();jthvec.clear();jphvec.clear();jvhvec.clear();
			//	getRecoverData(EDGECONTRACT, jpointvec, jthvec, jphvec, jvhvec);
			//	optjournal->pushJournal(EDGECONTRACT, jpointvec, jthvec, jphvec, jvhvec);

			//	/* calculate the worst quality after edge contraction */
			//	pointstar.clear();
			//	tmesh->point_star(newpoint, pointstar);
			//	minqualnew = qualcalculator->minstackquality(tmesh, pointstar, qualityMetric);

			//	/* if the quality gets worse, undo */
			//	if (!(minqualnew > minqualbefore))
			//	{
			//		contracted = false;
			//		/* set recover data */
			//		record = optjournal->topJournal();
			//		optjournal->popJournal();
			//		setRecoverData(EDGECONTRACT, record.pointvec, record.thvec, record.phvec, record.vhvec);
			//		/* do recover */
			//		edgestar.clear();
 		//			recover_edge_contract(&edgestar);
			//		minqualnew = qualcalculator->minstackquality(tmesh, edgestar, qualityMetric);

			//		/* repair the invalid tetrahandle */
			//		tetrahandle = edgestar[0];
			//	}
			//}

			///* return with success */
			//if (contracted) return true;

			///* if we're requested just to try the first edge, bail now */
			//if (justfirstedge) return false;
		}
		return false;
	}

	/************************************************************************/
	/* Recover data read and write                                          */
	/************************************************************************/

	/* get a piece of journal record */
	void TMeshToperator::getRecoverData(int type, std::vector<Point> &pointvec, std::vector<TetraHandle> &thvec, 
		                                std::vector<PointHandle> &phvec, std::vector<VertexHandle> &vhvec)
	{
		/* type specific stuff */
		switch (type)
		{
		case FLIP23:
			phvec.push_back(rec_flip23_top);
			phvec.push_back(rec_flip23_bot);
			phvec.push_back(rec_flip23_p[0]);
			phvec.push_back(rec_flip23_p[1]);
			phvec.push_back(rec_flip23_p[2]);
			thvec.push_back(rec_flip23_old_tetra[0]);
			thvec.push_back(rec_flip23_old_tetra[1]);
			thvec.push_back(rec_flip23_tetra[0]);
			thvec.push_back(rec_flip23_tetra[1]);
			thvec.push_back(rec_flip23_tetra[2]);
			break;
		case FLIP22:
			phvec.push_back(rec_flip22_top);
			phvec.push_back(rec_flip22_p1);
			phvec.push_back(rec_flip22_p2);
			phvec.push_back(rec_flip22_p[0]);
			phvec.push_back(rec_flip22_p[1]);
			thvec.push_back(rec_flip22_old_tetra[0]);
			thvec.push_back(rec_flip22_old_tetra[1]);
			thvec.push_back(rec_flip22_tetra[0]);
			thvec.push_back(rec_flip22_tetra[1]);
			break;
		case FLIP32:
			phvec.push_back(rec_flip32_top);
			phvec.push_back(rec_flip32_bot);
			phvec.push_back(rec_flip32_p[0]);
			phvec.push_back(rec_flip32_p[1]);
			phvec.push_back(rec_flip32_p[2]);
			thvec.push_back(rec_flip32_old_tetra[0]);
			thvec.push_back(rec_flip32_old_tetra[1]);
			thvec.push_back(rec_flip32_old_tetra[2]);
			thvec.push_back(rec_flip32_tetra[0]);
			thvec.push_back(rec_flip32_tetra[1]);
			break;
		case FLIP13:
			phvec.push_back(rec_flip13_p[0]);
			phvec.push_back(rec_flip13_p[1]);
			phvec.push_back(rec_flip13_p[2]);
			phvec.push_back(rec_flip13_p[3]);
			thvec.push_back(rec_flip13_old_tetra);
			thvec.push_back(rec_flip13_tetra[0]);
			thvec.push_back(rec_flip13_tetra[1]);
			thvec.push_back(rec_flip13_tetra[2]);
			break;
		case FLIP12:
			phvec.push_back(rec_flip12_p[0]);
			phvec.push_back(rec_flip12_p[1]);
			phvec.push_back(rec_flip12_p[2]);
			phvec.push_back(rec_flip12_p[3]);
			thvec.push_back(rec_flip12_old_tetra);
			thvec.push_back(rec_flip12_tetra[0]);
			thvec.push_back(rec_flip12_tetra[1]);
			break;
		case FLIP14:
			phvec.push_back(rec_flip14_p[0]);
			phvec.push_back(rec_flip14_p[1]);
			phvec.push_back(rec_flip14_p[2]);
			phvec.push_back(rec_flip14_p[3]);
			thvec.push_back(rec_flip14_old_tetra);
			thvec.push_back(rec_flip14_tetra[0]);
			thvec.push_back(rec_flip14_tetra[1]);
			thvec.push_back(rec_flip14_tetra[2]);
			thvec.push_back(rec_flip14_tetra[3]);
			break;
		case FLIP41:
			pointvec.push_back(rec_flip41_midpoint);
			thvec.push_back(rec_flip41_old_tetra[0]);
			thvec.push_back(rec_flip41_old_tetra[1]);
			thvec.push_back(rec_flip41_old_tetra[2]);
			thvec.push_back(rec_flip41_old_tetra[3]);
			thvec.push_back(rec_flip41_tetra);
			break;
		case FLIP31:
			pointvec.push_back(rec_flip31_midpoint);
			phvec.push_back(rec_flip31_top);
			phvec.push_back(rec_flip31_face_p[0]);
			phvec.push_back(rec_flip31_face_p[1]);
			phvec.push_back(rec_flip31_face_p[2]);
			thvec.push_back(rec_flip31_old_tetra[0]);
			thvec.push_back(rec_flip31_old_tetra[1]);
			thvec.push_back(rec_flip31_old_tetra[2]);
			thvec.push_back(rec_flip31_tetra);
			break;
		case EDGECONTRACT:
			pointvec.push_back(rec_ec_erase_point);
			pointvec.push_back(rec_ec_point);
			phvec.push_back(rec_ec_point_handle);
			copy(rec_ec_tetravertex.begin(), rec_ec_tetravertex.end(), back_inserter(phvec));
			copy(rec_ec_old_tetra.begin(), rec_ec_old_tetra.end(), back_inserter(thvec));
			copy(rec_ec_to_point_vertex_container.begin(), rec_ec_to_point_vertex_container.end(), back_inserter(thvec));
			break;
		default:
			printf("Undefined operation type %d\n", type);
			exit(1);
		};

	}

	/* set operation recover data with a piece of journal record */
	void TMeshToperator::setRecoverData(int type, std::vector<Point> &pointvec, std::vector<TetraHandle> &thvec, 
		                                std::vector<PointHandle> &phvec, std::vector<VertexHandle> &vhvec)
	{
		/* type specific stuff */
		switch (type)
		{
		case FLIP23:
			if (phvec.size()==5&&thvec.size()==5)
			{
				rec_flip23_valid = true;
				rec_flip23_top = phvec[0];
				rec_flip23_bot = phvec[1];
				rec_flip23_p[0] = phvec[2];
				rec_flip23_p[1] = phvec[3];
				rec_flip23_p[2] = phvec[4];
				rec_flip23_old_tetra[0] = thvec[0];
				rec_flip23_old_tetra[1] = thvec[1];
				rec_flip23_tetra[0] = thvec[2];
				rec_flip23_tetra[1] = thvec[3];
				rec_flip23_tetra[2] = thvec[4];
			}
			else
				rec_flip23_valid = false;
			break;
		case FLIP22:
			if (phvec.size()==5 && thvec.size()==4)
			{
				rec_flip22_valid = true;
				rec_flip22_top = phvec[0];
				rec_flip22_p1 = phvec[1];
				rec_flip22_p2 = phvec[2];
				rec_flip22_p[0] = phvec[3];
				rec_flip22_p[1] = phvec[4];
				rec_flip22_old_tetra[0] = thvec[0];
				rec_flip22_old_tetra[1] = thvec[1];
				rec_flip22_tetra[0] = thvec[2];
				rec_flip22_tetra[1] = thvec[3];
			}
			else
				rec_flip22_valid = false;
			break;
		case FLIP32:
			if (phvec.size()==5 && thvec.size()==5)
			{
				rec_flip32_valid = true;
				rec_flip32_top = phvec[0];
				rec_flip32_bot = phvec[1];
				rec_flip32_p[0] = phvec[2];
				rec_flip32_p[1] = phvec[3];
				rec_flip32_p[2] = phvec[4];
				rec_flip32_old_tetra[0] = thvec[0];
				rec_flip32_old_tetra[1] = thvec[1];
				rec_flip32_old_tetra[2] = thvec[2];
				rec_flip32_tetra[0] = thvec[3];
				rec_flip32_tetra[1] = thvec[4];
			}
			else
				rec_flip32_valid = false;
			break;
		case FLIP13:
			if (phvec.size()==4 && thvec.size()==4)
			{
				rec_flip13_valid = true;
				rec_flip13_p[0] = phvec[0];
				rec_flip13_p[1] = phvec[1];
				rec_flip13_p[2] = phvec[2];
				rec_flip13_p[3] = phvec[3];
				rec_flip13_old_tetra = thvec[0];
				rec_flip13_tetra[0] = thvec[1];
				rec_flip13_tetra[1] = thvec[2];
				rec_flip13_tetra[2] = thvec[3];
			}
			else
				rec_flip13_valid = false;
			break;
		case FLIP12:
			if (phvec.size()==4 && thvec.size()==3)
			{
				rec_flip12_valid = true;
				rec_flip12_p[0] = phvec[0];
				rec_flip12_p[1] = phvec[1];
				rec_flip12_p[2] = phvec[2];
				rec_flip12_p[3] = phvec[3];
				rec_flip12_old_tetra = thvec[0];
				rec_flip12_tetra[0] = thvec[1];
				rec_flip12_tetra[1] = thvec[2];
			}
			else
				rec_flip12_valid = false;
			break;
		case FLIP14:
			if (phvec.size()==4 && thvec.size()==5)
			{
				rec_flip14_valid = true;
				rec_flip14_p[0] = phvec[0];
				rec_flip14_p[1] = phvec[1];
				rec_flip14_p[2] = phvec[2];
				rec_flip14_p[3] = phvec[3];
				rec_flip14_old_tetra = thvec[0];
				rec_flip14_tetra[0] = thvec[1];
				rec_flip14_tetra[1] = thvec[2];
				rec_flip14_tetra[2] = thvec[3];
				rec_flip14_tetra[3] = thvec[4];
			}
			else
				rec_flip14_valid = false;
			break;
		case FLIP41:
			if (pointvec.size()==1 && thvec.size()==5)
			{
				rec_flip41_valid = true;
				rec_flip41_midpoint = pointvec[0];
				rec_flip41_old_tetra[0] = thvec[0];
				rec_flip41_old_tetra[1] = thvec[1];
				rec_flip41_old_tetra[2] = thvec[2];
				rec_flip41_old_tetra[3] = thvec[3];
				rec_flip41_tetra = thvec[4];
			}
			else
				rec_flip41_valid = false;
			break;
		case FLIP31:
			if (pointvec.size()==1 && phvec.size()==4 && thvec.size()==4)
			{
				rec_flip31_valid = true;
				rec_flip31_midpoint = pointvec[0];
				rec_flip31_top = phvec[0];
				rec_flip31_face_p[0] = phvec[1];
				rec_flip31_face_p[1] = phvec[2];
				rec_flip31_face_p[2] = phvec[3];
				rec_flip31_old_tetra[0] = thvec[0];
				rec_flip31_old_tetra[1] = thvec[1];
				rec_flip31_old_tetra[2] = thvec[2];
				rec_flip31_tetra = thvec[3];
			}
			else
				rec_flip31_valid = false;
			break;
		case EDGECONTRACT:
			if (pointvec.size()==2)
			{
				int tetcnt;
				tetcnt = (phvec.size() - 1) / 2;
				rec_ec_valid = true;
				rec_ec_erase_point = pointvec[0];
				rec_ec_point = pointvec[1];
				rec_ec_point_handle = phvec[0];
				rec_ec_tetravertex.clear();
				copy(phvec.begin()+1, phvec.end(),back_inserter(rec_ec_tetravertex));
				rec_ec_old_tetra.clear();
				copy(thvec.begin(), thvec.begin()+tetcnt, back_inserter(rec_ec_old_tetra));
				rec_ec_to_point_vertex_container.clear();
				copy(thvec.begin()+tetcnt, thvec.end(), back_inserter(rec_ec_to_point_vertex_container));
			}
			else
				rec_ec_valid = false;
			break;
		default:
			printf("Undefined operation type %d\n", type);
			exit(1);
		};

	}

	/* operation recover */
	bool TMeshToperator::recover(int optType, std::vector<TetraHandle> *tetraVec_)
	{
		int flag;
		switch (optType)
		{
		case FLIP23:
			flag = recover_flip23(tetraVec_);
			break;
		case FLIP22:
			flag = recover_flip22(tetraVec_);
			break;
		case FLIP32:
			flag = recover_flip32(tetraVec_);
			break;
		case FLIP13:
			flag = recover_flip13(tetraVec_);
			break;
		case FLIP12:
			flag = recover_flip12(tetraVec_);
			break;
		case FLIP14:
			flag = recover_flip14(tetraVec_);
			break;
		case FLIP41:
			flag = recover_flip41(tetraVec_);
			break;
		case FLIP31:
			flag = recover_flip31(tetraVec_);
			break;
		case EDGECONTRACT:
			flag = recover_edge_contract(tetraVec_);
			break;
		default:
			printf("Undefined operation type %d\n", optType);
			exit(1);
		}
		if (flag)
			return false;
		return true;
	}

}