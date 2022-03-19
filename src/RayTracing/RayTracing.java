package RayTracing;

import javax.swing.*;
import java.awt.*;
import java.util.Objects;

class Constants {
    public static int canvasHeight = 600;
    public static int canvasWidth = 800;
    public static int viewPortHeight = 1;
    public static int viewPortWidth = 1;
}

class Vector3d {
    double x, y, z;

    public Vector3d() {
        x = y = z = 0;
    }

    public Vector3d(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public static double Dot(Vector3d v1, Vector3d v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    public static double Length(Vector3d v) {
        return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    public static Vector3d MultMatrixVect(Vector3d v, double[][] matrixRot) {
        return new Vector3d(matrixRot[0][0] * v.x + matrixRot[0][1] * v.y + matrixRot[0][2] * v.z,
                matrixRot[1][0] * v.x + matrixRot[1][1] * v.y + matrixRot[1][2] * v.z,
                matrixRot[2][0] * v.x + matrixRot[2][1] * v.y + matrixRot[2][2] * v.z);
    }
}

class Vector2d {
    double x, y;

    public Vector2d() {
        x = y = 0;
    }

    public Vector2d(double x, double y) {
        this.x = x;
        this.y = y;
    }
}

class Sphere_t {
    Sphere sphere;
    double t;

    public Sphere_t(Sphere sphere, double t) {
        this.sphere = sphere;
        this.t = t;
    }
}

class Sphere {
    Vector3d center;
    Color color;
    double radius;
    int specular; //зеркальность

    public Sphere(Vector3d center, double radius, Color color, int specular) {
        this.center = center;
        this.radius = radius;
        this.color = color;
        this.specular = specular;
    }
}

class Plane {
    Vector3d point;
    Vector3d n;
    Color color;

    public Plane(Vector3d point, Vector3d n, Color color) {
        this.point = point;
        this.n = n;
        this.color = color;
    }
}

class Light {
    String type;
    float intensity;
    Vector3d position;
    Vector3d direction;

    public Light(String type, float intensity, Vector3d position, Vector3d direction) {
        this.type = type;
        this.intensity = intensity;
        this.position = position;
        this.direction = direction;
    }
}

class DrawComponent extends JComponent {
    public void paint(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        Vector3d O = new Vector3d(0, 0, 0); //положение камеры
        double projection_plane_d = 1;
        double degrees = 0; //уголь поворота
        double radians = Math.toRadians(degrees);
        //матрицы поворота
        double Mupdown[][] = {{1, 0, 0}, {0, Math.cos(radians), -Math.sin(radians)}, {0, Math.sin(radians), Math.cos(radians)}};
        double Mrightleft[][] = {{Math.cos(radians), 0, Math.sin(radians)}, {0, 1, 0}, {-Math.sin(radians), 0, Math.cos(radians)}};
        //инициализация наших объектов
        Sphere[] spheres = {
                new Sphere(new Vector3d(0, 0, 8), 1, new Color(255, 0, 0), 500),
                new Sphere(new Vector3d(2, 0, 10), 1, new Color(0, 0, 255), 500),
                new Sphere(new Vector3d(-2, 0, 4), 1, new Color(0, 255, 0), 500)
        };
        Plane[] planes = {
                new Plane(new Vector3d(0, 0, 5), new Vector3d(0, 0, 1), new Color(200, 162, 200)),
                //new Plane(new Vector3d(0,0,10), new Vector3d(0,0,1), new Color(200,162,200))
        };
        Light[] lights = {
                new Light("ambient", 0.2f, null, null),
                new Light("point", 0.6f, new Vector3d(0, 0, 4), null),
                new Light("directional", 0.2f, null, new Vector3d(1, 4, 4))
        };
        for (int S_x = 0; S_x <= Constants.canvasWidth; S_x++) {
            int x = S_x - Constants.canvasWidth / 2;
            for (int S_y = 0; S_y <= Constants.canvasHeight; S_y++) {
                int y = Constants.canvasHeight / 2 - S_y;
                Vector3d D = CanvasToViewport(x, y, projection_plane_d);
                Vector3d nD = Vector3d.MultMatrixVect(D, Mupdown);
                //для сфер
                Color color = TraceRay(O, nD, 1, Double.POSITIVE_INFINITY, spheres, lights);
                //для плоскости
                //Color color = TraceRay_Plane(O, nD, 1, Double.POSITIVE_INFINITY, planes, lights);
                g2.setPaint(color);
                g2.fillRect(S_x, S_y, 1, 1);
            }
        }
    }

    private Vector3d CanvasToViewport(double x, double y, double d) { //d-расстояние до плоскости проекции
        return new Vector3d(x * Constants.viewPortWidth / Constants.canvasWidth, y * Constants.viewPortHeight / Constants.canvasHeight, d);
    }

    //вычисляет пересечение луча с каждой сферой и возвращает цвет сферы в ближайшей точке пересечения
    private Color TraceRay(Vector3d O, Vector3d D, double tMin, double tMax, Sphere[] spheres, Light[] lights) {
        Sphere_t st = ClosestIntersection(O, D, tMin, tMax, spheres);
        Sphere closest_sphere = st.sphere;
        double closest_t = st.t;
        if (closest_sphere == null) {
            return getBackground();
        }
        Vector3d P = new Vector3d(O.x + closest_t * D.x, O.y + closest_t * D.y, O.z + closest_t * D.z); //вычисление пересечения
        //вычисление нормали сферы в точке пересечения
        Vector3d N = new Vector3d(P.x - closest_sphere.center.x, P.y - closest_sphere.center.y, P.z - closest_sphere.center.z);
        Vector3d N_real = new Vector3d(N.x / Vector3d.Length(N), N.y / Vector3d.Length(N), N.z / Vector3d.Length(N));
        int argb = closest_sphere.color.getRGB();
        int alpha = (argb >> 24) & 0xff;
        int red = (argb >> 16) & 0xff;
        int green = (argb >> 8) & 0xff;
        int blue = (argb) & 0xff;
        float r, g, b;
        r = (float) red / 255;
        g = (float) green / 255;
        b = (float) blue / 255;
        Vector3d _D = new Vector3d(-D.x, -D.y, -D.z);
        float f_red = r * ComputeLighting(P, N_real, _D, closest_sphere.specular, lights, spheres);
        float f_green = g * ComputeLighting(P, N_real, _D, closest_sphere.specular, lights, spheres);
        float f_blue = b * ComputeLighting(P, N_real, _D, closest_sphere.specular, lights, spheres);
        if (f_red > 1)
            f_red = 1;
        if (f_green > 1)
            f_green = 1;
        if (f_blue > 1)
            f_blue = 1;
        return new Color(f_red, f_green, f_blue);
    }

    private Sphere_t ClosestIntersection(Vector3d O, Vector3d D, double tMin, double tMax, Sphere[] spheres) {
        double closest_t = Double.POSITIVE_INFINITY;
        Sphere closest_sphere = null;
        for (Sphere sphere :
                spheres) {
            Vector2d t1andt2 = IntersectRaySphere(O, D, sphere);
            if (t1andt2.x >= tMin && t1andt2.x <= tMax && t1andt2.x < closest_t) {
                closest_t = t1andt2.x;
                closest_sphere = sphere;
            }
            if (t1andt2.y >= tMin && t1andt2.y <= tMax && t1andt2.y < closest_t) {
                closest_t = t1andt2.y;
                closest_sphere = sphere;
            }
        }
        return new Sphere_t(closest_sphere, closest_t);
    }

    private Vector2d IntersectRaySphere(Vector3d O, Vector3d D, Sphere s) {
        Vector3d C = s.center;
        double r = s.radius;
        Vector3d oc = new Vector3d(O.x - C.x, O.y - C.y, O.z - C.z);

        double k1, k2, k3, disc, t1, t2;
        k1 = Vector3d.Dot(D, D);
        k2 = 2 * Vector3d.Dot(oc, D);
        k3 = Vector3d.Dot(oc, oc) - r * r;

        disc = k2 * k2 - 4 * k1 * k3;
        if (disc < 0)
            return new Vector2d(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);

        t1 = (-k2 + Math.sqrt(disc)) / (2 * k1);
        t2 = (-k2 - Math.sqrt(disc)) / (2 * k1);
        return new Vector2d(t1, t2);
    }

    private float ComputeLighting(Vector3d P, Vector3d N, Vector3d V, int s, Light[] lights, Sphere[] spheres) {
        float i = 0;
        Vector3d L;
        double tMax;
        for (Light light : lights) {
            if (Objects.equals(light.type, "ambient")) {
                i += light.intensity;
            } else {
                if (Objects.equals(light.type, "point")) {
                    L = new Vector3d(light.position.x - P.x, light.position.y - P.y, light.position.z - P.z);
                    tMax = 1;
                } else {
                    L = new Vector3d(light.direction.x, light.direction.y, light.direction.z);
                    tMax = Double.POSITIVE_INFINITY;
                }
                double n_dot_l = Vector3d.Dot(N, L);
                //проверка тени
                Sphere_t st = ClosestIntersection(P, L, 0.001, tMax, spheres);
                if (st.sphere != null)
                    continue;
                //диффузность
                if (n_dot_l > 0)
                    i += light.intensity * n_dot_l / (Vector3d.Length(N) * Vector3d.Length(L));
                //зеркальность
                if (s != 1) {
                    Vector3d R = new Vector3d(2 * Vector3d.Dot(N, L) * N.x - L.x,
                            2 * Vector3d.Dot(N, L) * N.y - L.y,
                            2 * Vector3d.Dot(N, L) * N.z - L.z);
                    double r_dot_v = Vector3d.Dot(R, V);
                    if (r_dot_v > 0) {
                        i += light.intensity * Math.pow(r_dot_v / (Vector3d.Length(R) * Vector3d.Length(V)), s);
                    }
                }
            }
        }
        return i;
    }

}

class MyFrame extends JFrame {

    public MyFrame() {
        setSize(Constants.canvasWidth, Constants.canvasHeight);
        setTitle("MyPicture");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setVisible(true);
        DrawComponent dc = new DrawComponent();
        add(dc);
    }
}

public class RayTracing {
    public static void main(String[] args) {
        new MyFrame();
    }
}
