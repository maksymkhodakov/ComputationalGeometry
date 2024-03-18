module org.example.graphics {
    requires javafx.controls;
    requires javafx.fxml;
    requires jts;


    opens org.example.graphics to javafx.fxml;
    exports org.example.graphics;
}