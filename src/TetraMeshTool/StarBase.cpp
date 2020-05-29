#include <TetraMeshTool/StarBase.h>
#include <fstream>
#include <Windows.h>
namespace VolumeMesh
{
	/* check if a facet (indicated by fp1, fp2, fp3) is oriented to a vertex v */
	int orient(Point v, Point fp1, Point fp2, Point fp3)
	{
		Point center;
		TetraMesh::Normal facenormal;
		double result;
		center = (fp1 + fp2 + fp3) / 3.0;
		facenormal = (fp2 - fp1) % (fp3 - fp1);
		facenormal.normalize();
		result = (v - center) | facenormal;
		if (result > 0)
			return 1;
		if (result == 0)
			return 0;
		return -1;
	}

	/* get four points of a tetrahedron indexed by a TetraHandle*/
	bool tetPoints(TetraMesh *mesh, TetraHandle handle, Point *p)
	{
		if (p == NULL)
			p = new Point[4];
		if (handle.idx()<0 || handle.idx() >= (int)mesh->size_tetrahedron() )
			return false;
		TetraMesh::HedronVertexIter hv_iter;
		hv_iter = mesh->hedron_vertex_iter(handle);
		p[0] = mesh->point(hv_iter.handle()); ++ hv_iter;
		p[1] = mesh->point(hv_iter.handle()); ++ hv_iter;
		p[2] = mesh->point(hv_iter.handle()); ++ hv_iter;
		p[3] = mesh->point(hv_iter.handle());
		return true;
	}

	/* get three points of a halfface */
	bool facePoints(TetraMesh *mesh, HalfFaceHandle handle, Point *p)
	{
		if (p == NULL)
			p = new Point[3];
		if (handle.idx()<0 || handle.idx() >= (int)mesh->size_halfface())
			return false;
		TetraMesh::HalfFaceVertexIter fv_iter;
		fv_iter = mesh->half_face_vertex_iter(handle);
		p[0] = mesh->point(fv_iter.handle()); ++ fv_iter;
		p[1] = mesh->point(fv_iter.handle()); ++ fv_iter;
		p[2] = mesh->point(fv_iter.handle());
		return true;
	}

	/* get the center point of a face */
	Point faceCenter(TetraMesh *mesh, TetraMesh::HalfFaceHandle fh)
	{
		Point center;
		Point p[3];
		TetraMesh::HalfFaceVertexIter fv_iter;
		fv_iter = mesh->half_face_vertex_iter(fh);
		p[0] = mesh->point(fv_iter.handle()); ++fv_iter;
		p[1] = mesh->point(fv_iter.handle()); ++fv_iter;
		p[2] = mesh->point(fv_iter.handle());
		center = (p[0] + p[1] + p[2])/3.0;
		return center;
	}


	/* get tetra volume */
	double tetVolume(Point p1, Point p2, Point p3, Point p4)
	{
		return (((p2 - p1) % (p3 - p1)) | (p4 - p1)) / 6.0;
	}

	/* Point Shading (group) */
	void pointGrouping(int **neighbours, int *neighbourcnt, int pointcnt, int *pointcolour, int &colourcnt, float &time)
	{
		if (neighbours == NULL || neighbourcnt == NULL || pointcolour == NULL)
			return;

		LARGE_INTEGER li_start, li_finish;
		colourcnt = 0;
		int colourvec[50];
		int npidx;

		QueryPerformanceCounter(&li_start);
		// initialize the pointgroup array
		for (int i = 0; i < pointcnt; i++)
		{
			pointcolour[i] = -1;
		}

		// for each point
		for (int i = 0; i < pointcnt; i++)
		{
			memset(colourvec, 0, sizeof(int)*colourcnt);
			// record the color of neighbor points
			for (int j = 0; j < neighbourcnt[i]; j++)
			{
				npidx = neighbours[i][j];
				if (pointcolour[npidx] != -1)
				{
					colourvec[pointcolour[npidx]] = 1;
				}
			}
			// search the exist unused color
			for (int j = 0; j < colourcnt; j++)
			{
				if (!colourvec[j])
				{
					pointcolour[i] = j;
					break;
				}
			}
			// if exist color all used then add a new color
			if (pointcolour[i] == -1)
			{
				pointcolour[i] = colourcnt++;
			}
		}
		QueryPerformanceCounter(&li_finish);
		time = 0.0;
		time += float(li_finish.QuadPart - li_start.QuadPart);




		//// test
		//std::ofstream outfile;
		//outfile.open("F:\\testfile.txt", std::fstream::app);
		//int *pcarray;
		//pcarray = new int[colourcnt];
		//memset(pcarray,0,colourcnt*sizeof(int));
		//for (int i=0;i<pointcnt;i++)
		//{
		//	++ pcarray[pointcolour[i]];
		//}
		//outfile << colourcnt << "\n";
		//for (int i = 0; i < colourcnt; i++)
		//{
		//	outfile << "Colour " << i << " : " << pcarray[i] << "\n";
		//}
		//outfile << "\n";
		//outfile.close();
	}
}