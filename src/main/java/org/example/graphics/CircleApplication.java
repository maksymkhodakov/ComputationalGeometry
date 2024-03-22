package org.example.graphics;

import com.vividsolutions.jts.geom.*;
import com.vividsolutions.jts.triangulate.VoronoiDiagramBuilder;
import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
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
    private final Map<Long, Integer> mapStatistics = new HashMap<>();
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
        Button showStatsButton = new Button("Показати статистику");
        controls.getChildren().add(showStatsButton);
        showStatsButton.setOnAction(e -> showStatistics());
        controls.getChildren().addAll(new Label("n: "), nField, new Label("m: "), mField, new Label("points: "), numberOfPointsField, generateButton, radiusLabel);

        root.setTop(controls);
        setupCanvasZoom();
        setupDraggableCanvas(canvas);
        root.setCenter(canvas);

        generateButton.setOnAction(e -> {
            long startTime = System.currentTimeMillis();
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

            Polygon starPolygon = createStarShapePolygon(starShape); // starShape це список точок зіркового многокутника
            Geometry voronoiDiagram = generateVoronoiDiagram(starShape); // точки які використовуються для генерації діаграм
            drawVoronoiDiagram(voronoiDiagram, starPolygon); // Now pass the starPolygon as well
            final List<LineString> voronoiEdges = captureVoronoiEdges(voronoiDiagram);
            System.out.println("All Voronoi diagram's edges");
            System.out.println(voronoiEdges);

            final Set<Point> intersections = findAllIntersectionPoints(voronoiEdges);

            GraphicsContext gc = canvas.getGraphicsContext2D();
            gc.setFill(Color.GREEN);
            for (Point intersection : intersections) {
                gc.fillOval(intersection.x - 3, intersection.y - 3, 6, 6); // Малюємо маленький кружечок
            }
            System.out.println("Intersections:");
            System.out.println(intersections);


            // зберігаємо кружечки в списку
            List<Circle> circles = new ArrayList<>();
            for (Point intersection : intersections) {
                Polygon bufferedStarPolygon = (Polygon) starPolygon.buffer(0.01); // налаштовуємо буфферизацію для фільтрації точок
                if (bufferedStarPolygon.contains(new GeometryFactory().createPoint(new Coordinate(intersection.x, intersection.y)))) {
                    double radius = calculateDistanceToStarPolygonBoundary(intersection, starPolygon);
                    circles.add(new Circle(intersection, radius));
                }
            }

            // малюємо кружечки
            gc.setStroke(Color.BLUE);
            for (Circle circle : circles) {
                double x = circle.center.x;
                double y = circle.center.y;
                double radius = circle.radius;
                gc.strokeOval(x - radius, y - radius, 2 * radius, 2 * radius);
            }

            // список всіх кружечків та радіусів
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
            long endTime = System.currentTimeMillis();
            mapStatistics.put(endTime - startTime, n);
        });


        stage.setScene(scene);
        stage.setTitle("Зірковий многогранник і Вписане Коло");
        stage.show();
    }

    private void showStatistics() {
        Stage stage = new Stage();

        // Create new axes and chart instances for the new window
        NumberAxis xAxis = new NumberAxis();
        NumberAxis yAxis = new NumberAxis();
        LineChart<Number, Number> lineChart = new LineChart<>(xAxis, yAxis);

        xAxis.setLabel("Number of Vertices");
        yAxis.setLabel("Execution Time (ms)");
        lineChart.setTitle("Time Dependency on Number of Vertices");

        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Execution Time");

        // Populate the series with data from mapStatistics
        mapStatistics.forEach((time, vertices) -> series.getData().add(new XYChart.Data<>(vertices, time)));

        lineChart.getData().add(series);

        Scene scene = new Scene(lineChart, 800, 600);
        stage.setScene(scene);
        stage.setTitle("Performance Statistics");
        stage.show();
    }

    private double calculateDistanceToStarPolygonBoundary(Point point, Polygon starPolygon) {
        GeometryFactory geometryFactory = new GeometryFactory();
        Geometry pointGeometry = geometryFactory.createPoint(new Coordinate(point.x, point.y));
        Geometry boundary = starPolygon.getExteriorRing();
        return pointGeometry.distance(boundary);
    }

    public Set<Point> findAllIntersectionPoints(List<LineString> lines) {
        Set<com.vividsolutions.jts.geom.Point> intersectionPoints = new HashSet<>();

        for (int i = 0; i < lines.size(); i++) {
            for (int j = i + 1; j < lines.size(); j++) {
                Geometry intersection = lines.get(i).intersection(lines.get(j));

                if (intersection != null && !intersection.isEmpty()) {
                    if (intersection instanceof com.vividsolutions.jts.geom.Point) {
                        intersectionPoints.add((com.vividsolutions.jts.geom.Point) intersection);
                    } else if (intersection instanceof MultiPoint) {
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
        gc.setStroke(Color.ORANGE);
        gc.setLineWidth(1);

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

        gc.setFill(Color.BLUE);
        for (Point point : points) {
            Geometry pointGeometry = new GeometryFactory().createPoint(new Coordinate(point.x, point.y));
            if (starPolygon.contains(pointGeometry)) {
                gc.fillOval(point.x - 3, point.y - 3, 6, 6);
            }
        }
    }

    private Polygon createStarShapePolygon(List<Point> starShapePoints) {
        GeometryFactory geometryFactory = new GeometryFactory();
        Coordinate[] coordinates = new Coordinate[starShapePoints.size() + 1];
        for (int i = 0; i < starShapePoints.size(); i++) {
            coordinates[i] = new Coordinate(starShapePoints.get(i).x, starShapePoints.get(i).y);
        }
        coordinates[starShapePoints.size()] = coordinates[0];
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
