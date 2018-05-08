using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp1
{
	class Program
	{
		static void Main(string[] args)
		{
			Simplify test = new Simplify();
			test.load_obj("bunny.obj");
			test.simplify_mesh(1000);
			test.write_obj("bunny2.obj");
			vec3f v = new vec3f();
			test.GetVertex(0, v);
			test.GetVertex(1, v);
			IntVector i = new IntVector(3);
			test.GetTriangle(0, i);
		}
	}
}
