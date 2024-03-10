package org.example.graphics;

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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

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
        controls.getChildren().addAll(numberOfPointsField, generateButton, radiusLabel);

        root.setTop(controls);
        setupCanvasZoom();
        setupDraggableCanvas(canvas);
        root.setCenter(canvas);

        generateButton.setOnAction(e -> {
            int numberOfPoints;
            try {
                numberOfPoints = Integer.parseInt(numberOfPointsField.getText());
            } catch (NumberFormatException ex) {
                // рандомно згенероване від 100 до 1000
                numberOfPoints = ThreadLocalRandom.current().nextInt(100, 1000);
            }
            generatePoints(numberOfPoints);
            drawPoints();
            List<Point> convexHull = findConvexHull(points);
            drawConvexHull(convexHull);
            Circle inscribedCircle = findInscribedCircle(convexHull);
            drawCircle(inscribedCircle);
            radiusLabel.setText(String.format("Радіус кола: %.2f", inscribedCircle.radius)); // Виведення радіуса на мітці
        });

        stage.setScene(scene);
        stage.setTitle("Опукла Оболонка і Вписане Коло");
        stage.show();
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

    private void drawCircle(Circle circle) {
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.setStroke(Color.BLUE);
        gc.setLineWidth(1);
        gc.strokeOval(circle.center.x - circle.radius, circle.center.y - circle.radius, 2 * circle.radius, 2 * circle.radius);
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

    private List<Point> findConvexHull(List<Point> points) {
        if (points.size() < 3) return Collections.emptyList();

        List<Point> sortedPoints = new ArrayList<>(points);
        sortedPoints.sort((p1, p2) -> p1.y != p2.y ? Double.compare(p1.y, p2.y) : Double.compare(p1.x, p2.x));
        Point p0 = sortedPoints.get(0);
        sortedPoints.sort((p1, p2) -> {
            double angle1 = Math.atan2(p1.y - p0.y, p1.x - p0.x);
            double angle2 = Math.atan2(p2.y - p0.y, p2.x - p0.x);
            return Double.compare(angle1, angle2);
        });

        List<Point> hull = new ArrayList<>();
        for (Point p : sortedPoints) {
            while (hull.size() >= 2 && crossProduct(hull.get(hull.size() - 2), hull.get(hull.size() - 1), p) <= 0) {
                hull.remove(hull.size() - 1);
            }
            hull.add(p);
        }

        return hull;
    }

    private double crossProduct(Point a, Point b, Point c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
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

    public Circle findInscribedCircle(List<Point> hull) {
        Point centroid = findCentroid(hull);
        double minRadius = Double.MAX_VALUE;

        for (int i = 0; i < hull.size(); i++) {
            Point start = hull.get(i);
            Point end = hull.get((i + 1) % hull.size());
            Line edge = new Line(start, end);
            double radius = distanceFromPointToLine(centroid, edge);
            if (radius < minRadius) {
                minRadius = radius;
            }
        }

        return new Circle(centroid, minRadius);
    }

    // Функція для обчислення відстані від точки до прямої
    private double distanceFromPointToLine(Point p, Line line) {
        return Math.abs(line.a * p.x + line.b * p.y + line.c) / Math.sqrt(line.a * line.a + line.b * line.b);
    }

    private Point findCentroid(List<Point> hull) {
        double centroidX = 0;
        double centroidY = 0;
        for (Point point : hull) {
            centroidX += point.x;
            centroidY += point.y;
        }
        return new Point(centroidX / hull.size(), centroidY / hull.size());
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

        Line(Point p1, Point p2) {
            a = p2.y - p1.y;
            b = p1.x - p2.x;
            c = -a * p1.x - b * p1.y;
        }
    }

    public static class Point {
        double x;
        double y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }
    }
}
