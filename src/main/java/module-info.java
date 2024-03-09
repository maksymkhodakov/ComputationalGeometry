module org.example.graphics {
    requires javafx.controls;
    requires javafx.fxml;


    opens org.example.graphics to javafx.fxml;
    exports org.example.graphics;
}