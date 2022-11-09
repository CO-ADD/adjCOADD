$(document).ready(function () {

    $(document).on("dblclick", ".editable", function () {
        console.log("loading editableTable...")

        var value = $(this).text();
        var input = "<input type='text' class='input-data' value='" + value + "' class='form-control'>";
        input += "<button class='btnclicksave' class='form-control'>save</button>";
        input += "<button class='btnclickcancel' class='form-control'>cancel</button>";
        $(this).html(input);
        $(this).removeClass("editable")
    });

    $(document).on("click", ".btnclicksave", function () {
        var value = $('.input-data').val();
        var td = $(this).parent("td");
        // $this.remove();
        td.html(value);
        td.addClass("editable");
        var type = td.data("type");
        sendToServer(td.data("id"), value, type);
    });
    $(document).on("click", ".btnclickcancel", function () {
        var value = $('.input-data').val();
        var td = $(this).parent("td");
        // $this.remove();
        td.html(value);
        td.addClass("editable");
    });



})
