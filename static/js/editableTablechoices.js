$(document).ready(function () {
    $(document).on("dblclick", ".editablechoices", function () {
        console.log("loading editableTable...")
        var selectbox = $(this).find(".select")
        var detail = $(this).find(".detail")
        selectbox.addClass("visible")
        $('.updateformshow').addClass("visible")
        detail.addClass("not-visible")
        $(this).removeClass("editablechoices")
    });

    $(document).keypress(
        function (event) {
            if (event.which == '13') {
                event.preventDefault();
            }
        });
})
