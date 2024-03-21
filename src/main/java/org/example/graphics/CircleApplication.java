package org.example.graphics;

import com.vividsolutions.jts.geom.*;
import com.vividsolutions.jts.triangulate.VoronoiDiagramBuilder;
import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.paint.Color;
import javafx.stage.Stage;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;

public class CircleApplication extends Application {
    private final List<Point> points = new ArrayList<>();
    private final Canvas canvas = new Canvas(1920, 1080);
    private final Label radiusLabel = new Label("Радіус кола: ");
    private double scaleFactor = 1.0; // Початковий масштаб
    private double lastX; // Координати для збереження попередньої позиції курсору
    private double lastY; // Координати для збереження попередньої позиції курсору


    @Override
    public void start(Stage stage) {
        BorderPane root = new BorderPane();
        Scene scene = new Scene(root);


        HBox controls = new HBox();
        TextField numberOfPointsField = new TextField();
        Button generateButton = new Button("Генерувати");
        TextField nField = new TextField(); // Для вводу кількості вершин
        TextField mField = new TextField(); // Для вводу кроку
        controls.getChildren().addAll(new Label("n: "), nField, new Label("m: "), mField, new Label("points: "), numberOfPointsField, generateButton, radiusLabel);

        root.setTop(controls);
        setupCanvasZoom();
        setupDraggableCanvas(canvas);
        root.setCenter(canvas);

        generateButton.setOnAction(e -> {
            int numberOfPoints;
            int n;
            int m;
            try {
                numberOfPoints = Integer.parseInt(numberOfPointsField.getText());
                n = Integer.parseInt(nField.getText());
                m = Integer.parseInt(mField.getText());
            } catch (NumberFormatException ex) {
                // Якщо введення некоректне, використовуйте стандартні значення
                numberOfPoints = ThreadLocalRandom.current().nextInt(100, 10000);
                n = 5; // Стандартне значення для n
                m = 2; // Стандартне значення для m
            }
            generatePoints(numberOfPoints);
            drawPoints();
            // Модифікація тут для використання n і m
            List<Point> starShape = formStarShape(points, n, m);
            drawConvexHull(starShape); // Малювання зіркового многокутника

            Polygon starPolygon = createStarShapePolygon(starShape); // starShape is your list of star shape points
            Geometry voronoiDiagram = generateVoronoiDiagram(starShape); // points are used to generate the diagram
            drawVoronoiDiagram(voronoiDiagram, starPolygon); // Now pass the starPolygon as well
            final List<LineString> voronoiEdges = captureVoronoiEdges(voronoiDiagram);
            System.out.println("All Voronoi diagram's edges");
            System.out.println(voronoiEdges);

            final Set<Point> intersections = findAllIntersectionPoints(voronoiEdges);

            GraphicsContext gc = canvas.getGraphicsContext2D();
            // Highlighting intersection points
            gc.setFill(Color.GREEN); // Choose a color for the intersection points
            for (Point intersection : intersections) {
                gc.fillOval(intersection.x - 3, intersection.y - 3, 6, 6); // Draw a small circle for each intersection point
            }
            System.out.println("Intersections:");
            System.out.println(intersections);


            // Store circles and their radiuses in a list
            List<Circle> circles = new ArrayList<>();
            for (Point intersection : intersections) {
                Polygon bufferedStarPolygon = (Polygon) starPolygon.buffer(0.01); // Adjust the buffer distance as needed
                if (bufferedStarPolygon.contains(new GeometryFactory().createPoint(new Coordinate(intersection.x, intersection.y)))) {
                    double radius = calculateDistanceToStarPolygonBoundary(intersection, starPolygon);
                    circles.add(new Circle(intersection, radius));
                }
            }

            // Draw circles with radiuses up to the star polygon
            gc.setStroke(Color.BLUE); // Set the color for the circles
            for (Circle circle : circles) {
                double x = circle.center.x;
                double y = circle.center.y;
                double radius = circle.radius;
                gc.strokeOval(x - radius, y - radius, 2 * radius, 2 * radius);
            }

            // List all circles and their radiuses
            System.out.println("All circles inside and their radiuses:");
            for (int i = 0; i < circles.size(); i++) {
                Circle circle = circles.get(i);
                System.out.println("Circle " + (i + 1) + ": Center(" + circle.center.x + ", " + circle.center.y + "), Radius: " + circle.radius);
            }

            circles.stream().max(Comparator.comparing(circle -> circle.radius))
                    .ifPresent(circle -> {
                        gc.setStroke(Color.RED);
                        double x = circle.center.x;
                        double y = circle.center.y;
                        double radius = circle.radius;
                        gc.strokeOval(x - radius, y - radius, 2 * radius, 2 * radius);
                        radiusLabel.setText(String.format("Радіус кола: %.2f", circle.radius));
                    });
        });


        stage.setScene(scene);
        stage.setTitle("Опукла Оболонка і Вписане Коло");
        stage.show();
    }

    // Calculate the distance from an intersection point to the closest point on the star polygon boundary
    private double calculateDistanceToStarPolygonBoundary(Point point, Polygon starPolygon) {
        GeometryFactory geometryFactory = new GeometryFactory();
        Geometry pointGeometry = geometryFactory.createPoint(new Coordinate(point.x, point.y));
        Geometry boundary = starPolygon.getExteriorRing();
        return pointGeometry.distance(boundary);
    }

    public Set<Point> findAllIntersectionPoints(List<LineString> lines) {
        Set<com.vividsolutions.jts.geom.Point> intersectionPoints = new HashSet<>();

        // Compare each line with every other line exactly once
        for (int i = 0; i < lines.size(); i++) {
            for (int j = i + 1; j < lines.size(); j++) {
                Geometry intersection = lines.get(i).intersection(lines.get(j));

                // If there's an intersection, process it
                if (intersection != null && !intersection.isEmpty()) {
                    if (intersection instanceof com.vividsolutions.jts.geom.Point) {
                        // If the intersection is a single point, add it directly
                        intersectionPoints.add((com.vividsolutions.jts.geom.Point) intersection);
                    } else if (intersection instanceof MultiPoint) {
                        // If the intersection is multiple points, add them all
                        for (int pointIndex = 0; pointIndex < intersection.getNumGeometries(); pointIndex++) {
                            intersectionPoints.add((com.vividsolutions.jts.geom.Point) intersection.getGeometryN(pointIndex));
                        }
                    }
                }
            }
        }

        return intersectionPoints.stream().map(this::convertJtsPointToCustomPoint).collect(Collectors.toSet());
    }

    private Point convertJtsPointToCustomPoint(com.vividsolutions.jts.geom.Point jtsPoint) {
        return new Point(jtsPoint.getCoordinate().x, jtsPoint.getCoordinate().y);
    }


    public List<LineString> captureVoronoiEdges(Geometry voronoiDiagram) {
        List<LineString> edges = new ArrayList<>();
        GeometryFactory geometryFactory = new GeometryFactory();

        for (int i = 0; i < voronoiDiagram.getNumGeometries(); i++) {
            Geometry geometry = voronoiDiagram.getGeometryN(i);
            if (geometry instanceof Polygon) {
                Polygon polygon = (Polygon) geometry;
                Coordinate[] coordinates = polygon.getExteriorRing().getCoordinates();
                for (int j = 0; j < coordinates.length - 1; j++) {
                    Coordinate start = coordinates[j];
                    Coordinate end = coordinates[j + 1];
                    edges.add(geometryFactory.createLineString(new Coordinate[]{start, end}));
                }
            }
        }

        return edges;
    }

    private void drawVoronoiDiagram(Geometry voronoiDiagram, Polygon starPolygon) {
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.setStroke(Color.ORANGE); // Set a color that distinguishes Voronoi edges
        gc.setLineWidth(1); // Line width for Voronoi edges

        for (int i = 0; i < voronoiDiagram.getNumGeometries(); i++) {
            Geometry geometry = voronoiDiagram.getGeometryN(i);
            if (geometry instanceof Polygon) {
                Polygon polygon = (Polygon) geometry;
                Coordinate[] coordinates = polygon.getExteriorRing().getCoordinates();
                for (int j = 0; j < coordinates.length - 1; j++) {
                    Coordinate start = coordinates[j];
                    Coordinate end = coordinates[j + 1];
                    gc.strokeLine(start.x, start.y, end.x, end.y);
                }
            }
        }

        // Now, let's mark all points that are inside the star polygon
        gc.setFill(Color.BLUE); // Set a color for points inside the star polygon
        for (Point point : points) { // Assuming 'points' are the Voronoi vertices or any points you want to check
            Geometry pointGeometry = new GeometryFactory().createPoint(new Coordinate(point.x, point.y));
            if (starPolygon.contains(pointGeometry)) {
                gc.fillOval(point.x - 3, point.y - 3, 6, 6); // Draw a small circle for each point inside
            }
        }
    }

    private Polygon createStarShapePolygon(List<Point> starShapePoints) {
        GeometryFactory geometryFactory = new GeometryFactory();
        Coordinate[] coordinates = new Coordinate[starShapePoints.size() + 1];
        for (int i = 0; i < starShapePoints.size(); i++) {
            coordinates[i] = new Coordinate(starShapePoints.get(i).x, starShapePoints.get(i).y);
        }
        coordinates[starShapePoints.size()] = coordinates[0]; // Close the polygon by repeating the first point at the end
        LinearRing linearRing = geometryFactory.createLinearRing(coordinates);
        return geometryFactory.createPolygon(linearRing, null);
    }


    private Geometry generateVoronoiDiagram(List<Point> points) {
        GeometryFactory geometryFactory = new GeometryFactory();
        GeometryCollection geometryCollection = convertPointsToGeometryCollection(points, geometryFactory);

        VoronoiDiagramBuilder voronoiBuilder = new VoronoiDiagramBuilder();
        voronoiBuilder.setSites(geometryCollection);
        return voronoiBuilder.getDiagram(geometryFactory);
    }

    private GeometryCollection convertPointsToGeometryCollection(List<Point> points, GeometryFactory geometryFactory) {
        Geometry[] geometries = new Geometry[points.size()];
        for (int i = 0; i < points.size(); i++) {
            geometries[i] = geometryFactory.createPoint(new Coordinate(points.get(i).x, points.get(i).y));
        }
        return new GeometryCollection(geometries, geometryFactory);
    }

    private void setupCanvasZoom() {
        canvas.setOnScroll(e -> {
            double delta = 1.2;
            double scale = e.getDeltaY() < 0 ? 1 / delta : delta;

            scaleFactor *= scale;
            canvas.setScaleX(scaleFactor);
            canvas.setScaleY(scaleFactor);

            e.consume();
        });
    }

    private void setupDraggableCanvas(Canvas canvas) {
        canvas.addEventFilter(MouseEvent.MOUSE_PRESSED, event -> {
            lastX = event.getSceneX();
            lastY = event.getSceneY();
        });

        canvas.addEventFilter(MouseEvent.MOUSE_DRAGGED, event -> {
            double offsetX = event.getSceneX() - lastX;
            double offsetY = event.getSceneY() - lastY;

            canvas.setTranslateX(canvas.getTranslateX() + offsetX);
            canvas.setTranslateY(canvas.getTranslateY() + offsetY);

            lastX = event.getSceneX();
            lastY = event.getSceneY();
        });
    }

    private void generatePoints(int numberOfPoints) {
        points.clear();
        for (int i = 0; i < numberOfPoints; i++) {
            double x = Math.random() * canvas.getWidth();
            double y = Math.random() * canvas.getHeight();
            points.add(new Point(x, y));
        }
    }

    private void drawPoints() {
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.clearRect(0, 0, canvas.getWidth(), canvas.getHeight());
        gc.setFill(Color.BLACK);
        for (Point point : points) {
            gc.fillOval(point.x - 2, point.y - 2, 4, 4);
        }
    }


    private List<Point> formStarShape(List<Point> points, int n, int m) {
        if (points.size() < n) return Collections.emptyList(); // Переконатися, що є достатньо точок

        // Знайти центральну точку для визначення початкової вершини
        Point center = findCentroid(points);

        // Сортування точок за відстанню від центру
        points.sort(Comparator.comparingDouble(p -> Math.hypot(p.x - center.x, p.y - center.y)));

        // Вибір n найближчих точок до центру як потенційних вершин зірки
        List<Point> selectedPoints = points.subList(0, n);

        // Сортування вибраних точок за кутом відносно центру
        selectedPoints.sort((p1, p2) -> {
            double angle1 = Math.atan2(p1.y - center.y, p1.x - center.x);
            double angle2 = Math.atan2(p2.y - center.y, p2.x - center.x);
            return Double.compare(angle1, angle2);
        });

        // Формування зіркового многокутника через вибір точок з кроком m
        List<Point> starShape = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            starShape.add(selectedPoints.get((i * m) % n));
        }

        return starShape;
    }

    private void drawConvexHull(List<Point> hull) {
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.setStroke(Color.RED);
        gc.setLineWidth(2);

        for (int i = 0; i < hull.size(); i++) {
            Point p1 = hull.get(i);
            Point p2 = hull.get((i + 1) % hull.size());
            gc.strokeLine(p1.x, p1.y, p2.x, p2.y);
        }
    }

    public Point findCentroid(List<Point> hull) {
        double totalArea = 0;
        double centroidX = 0;
        double centroidY = 0;

        for (int i = 1; i < hull.size() - 1; i++) {
            double area = triangleArea(hull.get(0), hull.get(i), hull.get(i + 1));
            totalArea += area;
            centroidX += (hull.get(0).x + hull.get(i).x + hull.get(i + 1).x) * area;
            centroidY += (hull.get(0).y + hull.get(i).y + hull.get(i + 1).y) * area;
        }

        if (totalArea == 0) return new Point(0, 0);

        centroidX /= (3.0 * totalArea);
        centroidY /= (3.0 * totalArea);
        return new Point(centroidX, centroidY);
    }

    private double triangleArea(Point a, Point b, Point c) {
        return 0.5 * Math.abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    }

    public static class Circle {
        Point center;
        double radius;

        Circle(Point center, double radius) {
            this.center = center;
            this.radius = radius;
        }
    }

    public static class Line {
        double a;
        double b;
        double c;

        public Line(double a, double b, double c) {
            this.a = a;
            this.b = b;
            this.c = c;
        }
    }

    public static class Point {
        double x;
        double y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public int hashCode() {
            long temp = Double.doubleToLongBits(x);
            int result = (int) (temp ^ (temp >>> 32));
            temp = Double.doubleToLongBits(y);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }

        @Override
        public String toString() {
            return "Point{" + "x=" + x + ", y=" + y + '}';
        }
    }
}
