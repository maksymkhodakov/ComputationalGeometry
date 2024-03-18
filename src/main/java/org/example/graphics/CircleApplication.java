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
import java.util.Comparator;
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
            Circle inscribedCircle = findInscribedCircle(starShape);
            drawCircle(inscribedCircle);
            radiusLabel.setText(String.format("Радіус кола: %.2f", inscribedCircle.radius));
            List<Line> bisectors = generateBisectors(starShape);
            drawBisectors(bisectors);
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

    private List<Line> generateBisectors(List<Point> convexHull) {
        List<Line> bisectors = new ArrayList<>();
        int n = convexHull.size();
        for (int i = 0; i < n; i++) {
            Point p0 = convexHull.get((i - 1 + n) % n); // Попередня точка
            Point p1 = convexHull.get(i);                // Поточна точка
            Point p2 = convexHull.get((i + 1) % n);      // Наступна точка

            // Вектор p0->p1
            double dx1 = p1.x - p0.x;
            double dy1 = p1.y - p0.y;

            // Вектор p1->p2
            double dx2 = p2.x - p1.x;
            double dy2 = p2.y - p1.y;

            // Нормалізація векторів
            double len1 = Math.sqrt(dx1 * dx1 + dy1 * dy1);
            double len2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
            dx1 /= len1;
            dy1 /= len1;
            dx2 /= len2;
            dy2 /= len2;

            // Середній вектор (напрямок бісектриси)
            double bx = dx1 + dx2;
            double by = dy1 + dy2;
            double blen = Math.sqrt(bx * bx + by * by);
            bx /= blen; // Нормалізація бісектриси
            by /= blen;

            double c = -(bx * p1.x + by * p1.y);
            bisectors.add(new Line(bx, by, c));
        }

        return bisectors;
    }

    private void drawBisectors(List<Line> bisectors) {
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.setStroke(Color.GREEN);
        gc.setLineWidth(1);

        for (Line bisector : bisectors) {
            double x1;
            double y1;
            double x2;
            double y2;

            if (bisector.b != 0) {
                x1 = 0;
                y1 = -bisector.c / bisector.b;
                x2 = canvas.getWidth();
                y2 = (-bisector.a * x2 - bisector.c) / bisector.b;
            } else {
                x1 = -bisector.c / bisector.a;
                y1 = 0;
                x2 = x1;
                y2 = canvas.getHeight();
            }

            gc.strokeLine(x1, y1, x2, y2);
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

    public Circle findInscribedCircle(List<Point> starShape) {
        if (starShape.isEmpty()) return null;

        // Початкове припущення для центру кола
        Point center = findStarCenter(starShape);

        double minDistanceToEdge = Double.MAX_VALUE;
        // Шукаємо мінімальну відстань від центру до країв
        for (int i = 0; i < starShape.size(); i++) {
            Point start = starShape.get(i);
            Point end = starShape.get((i + 1) % starShape.size());
            Line edge = new Line(start, end);

            double distanceToEdge = distanceFromPointToLine(center, edge);
            minDistanceToEdge = Math.min(minDistanceToEdge, distanceToEdge);
        }

        return new Circle(center, minDistanceToEdge);
    }

    private Point findStarCenter(List<Point> starShape) {
        // Спробуємо знайти більш точний центр для зірки
        double minX = Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE;
        double maxY = Double.MIN_VALUE;
        for (Point p : starShape) {
            minX = Math.min(minX, p.x);
            maxX = Math.max(maxX, p.x);
            minY = Math.min(minY, p.y);
            maxY = Math.max(maxY, p.y);
        }
        // Повертаємо точку посередині між мінімальною та максимальною координатами
        return new Point((minX + maxX) / 2, (minY + maxY) / 2);
    }


    // Функція для обчислення відстані від точки до прямої
    private double distanceFromPointToLine(Point p, Line line) {
        return Math.abs(line.a * p.x + line.b * p.y + line.c) / Math.sqrt(line.a * line.a + line.b * line.b);
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

        Line(Point p1, Point p2) {
            a = p2.y - p1.y;
            b = p1.x - p2.x;
            c = -a * p1.x - b * p1.y;
        }

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
    }
}
