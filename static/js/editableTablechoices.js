$(document).ready(function () {

    console.log("loading choiceseditable....")


    $(document).on("dblclick", ".editablechoices", function () {
        console.log("loading editableTable...")
        var selectbox = $(this).find("select")
        selectbox.addClass("visible")
        input = "<button class='btncancel' class='form-control'><i class='bi bi-x'></i></button>";
        input += "<button class='btnsave' class='form-control'><i class='bi bi-check'></i></button>";
        // input += selectbox
        if ($(this).find("button").length === 0) {
            console.log($(this).find("button"))

            $(this).append(input);
        }
        $(this).removeClass("editablechoices")
    });

    $("select").on("click", function () {
        if (!$(this).hasClass('selected')) {
            $(this).addClass('selected')
        }

    });

    $(document).on("click", ".btncancel", function () {
        var td = $(this).parent("td");
        // var value = td.val()
        td.find("select").removeClass("visible")

        $(this).remove();
        $(".btnsave").remove();
        // td.html(value)
        td.addClass("editablechoices");
        // window.location.reload()
    });

    $(document).on("click", ".btnsave", function () {
        var td = $(this).parent("td");
        td.find("select").removeClass("visible")
        // $(".btncancel").remove();
        var value = td.find("select").val()
        console.log(typeof ('value'))
        if (typeof (value) === 'object') {
            value = Object.values(value).toString()
        }
        sendToServer(td.data("id"), value, td.data("type"))
        $(this).remove();
        td.addClass("editablechoices");
        // window.location.reload()

    });



})
