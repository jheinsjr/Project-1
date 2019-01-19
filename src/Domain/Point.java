/*
 * 
 */
package Domain;

/**
 *
 * @author Joe Heins
 */
public class Point {

    //Instance Variables
    float[] point;

    public Point() {
        point = new float[2];

    }

    public Point(float x, float y) {
        point = new float[2];
        point[0] = x;
        point[1] = y;

    }

    public float getX() {
        return point[0];
    }

    public void setX(float x) {
        point[0] = x;
    }

    public float getY() {
        return point[1];
    }

    public void setY(float y) {
        point[1] = y;
    }

    public boolean equals(Point point) {
        if (point == this) {
            return true;
        }
        if (null == point) {
            return false;
        }
        if (point.getClass() != getClass()) {
            return false;
        }
        final Point other = (Point) point;
        if (Double.doubleToLongBits(this.getX())
                != Double.doubleToLongBits(other.getX())) {
            return false;
        }
        return Double.doubleToLongBits(this.getY())
                == Double.doubleToLongBits(other.getY());
    }

    public Point add(Point p) {
        Point refPoint = new Point();
        refPoint.setX(p.getX() + this.getX());
        refPoint.setY(p.getY() + this.getY());

        return refPoint;
    }

    /**
     *
     * @param p
     */
    public void combinePoints(Point p) {
        this.point[0] += p.getX();
        this.point[1] += p.getY();

    }

    public Point[] mulitiply(Point p) {
        Point[] refPoint = new Point[2];

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                refPoint[i].point[j] = this.point[i] * p.point[j];
            }
        }

        return refPoint;
    }

    public Point divide(Point p) {
        Point quotient = new Point();
        quotient.setX(this.getX() / p.getX());
        quotient.setY(this.getY() / p.getY());

        return quotient;
    }

    public Point difference(Point p) {
        Point difference = new Point();
        difference.setX(this.getX() - p.getX());
        difference.setY(this.getY() - p.getY());

        return difference;
    }

}
