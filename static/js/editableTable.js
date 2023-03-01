$(document).ready(function () {
    $(document).on("dblclick", ".editable", function () {
        var selectbox = $(this).find(".select")
        var detail = $(this).find(".detail")
        selectbox.addClass("visible")
        $('.updateformshow').addClass("visible")
        detail.addClass("not-visible")
        $(this).removeClass("editable")
    });

    $(document).keypress(
        function (event) {
            if (event.which == '13') {
                event.preventDefault();
            }
        });
})
