$(document).ready(function () {

    console.log("loading choiceseditable....")
    var value = []
    $(document).on("dblclick", ".editablechoices", function () {
        console.log("loading editableTable...")
        var selectbox = $(this).find("select")
        selectbox.toggleClass("visible")

        input = "<button class='btncancel' class='form-control'>cancel</button>";
        input += "<button class='btnok' class='form-control'>ok</button>";
        $(this).append(input);
        $(this).removeClass("editable")
    });
    $(document).on("click", "option", function () {
        value.push($(this).val());
        console.log(value)
    });

    $(document).on("click", ".btncancel", function () {
        var td = $(this).parent("td");
        td.find("select").removeClass("visible")
        $(this).remove();
        $(".btnok").remove();
        value = []
        td.addClass("editable");

    });

    $(document).on("click", ".btnok", function () {
        var td = $(this).parent("td");
        if (value == "" && typeof (td.find("select").val()) === "string") {
            value = td.find("select").val()
        }
        td.find("select").removeClass("visible")
        if (value == "") {

            alert("Strain_Type is null, press 'Cancel' discarding modifications")
        }
        $(this).remove();
        $(".btncancel").remove();
        input = "<button class='btnexit' class='form-control'>cancel</button>";
        input += "<button class='btnsave' class='form-control'>save</button>";
        td.text(value)
        td.append(input);
    });


    $(document).on("click", ".btnexit", function () {
        var td = $(this).parent("td");
        td.find("select").removeClass("visible")
        $(this).remove();
        $(".btnsave").remove();
        $(".btnok").remove();
        $(".btncancel").remove();
        value = []
        td.addClass("editable");
        window.location.reload()

    });

    $(document).on("click", ".btnsave", function () {
        var td = $(this).parent("td");
        td.find("select").removeClass("visible")
        $(this).remove();
        $(".btnexit").remove();
        $(".btnok").remove();
        $(".btncancel").remove();
        sendToServer(td.data("id"), value, td.data("type"))
        td.addClass("editable");

    });

    $("html").click(function () {
        //close popup
    });
})
