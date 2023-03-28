$(document).ready(function () {
  $(document).on("click", ".editablerow_stock", function () {
    console.log("editable stock...");
    var value = $(this).text();
    var input =
      "<p class='input-data'>Are you sure: " + value + " (n_left) - 1? </p>" + "<button data-value='" + value + "' class='btn_cancel btn btn-small border border-0'><i class='fa-solid fa-xmark'></i></button>" + "<button data-value='" + value + "' class='btn_ok btn btn-small border border-0'><i class='fa-solid fa-check'></i></button>";
    $(this).html(input);

    $(this).removeClass("editablerow_stock");
  });

  $(document).on("click", ".btn_ok", function () {
    console.log("button ok")
    var value = Number($(this).data('value')) - 1;
    var td = $(this).parent("td");
    $(this).remove();
    $(this).closest('button').remove();
    td.html(value.toString());
    td.addClass("editablerow_stock");
    var name = td.data("name");
    var type = td.data("type");
    data = { name: name, type: type, value: value };
    console.log(data);
    sendToServer(data, td);
  });
  $(document).on("click", ".btn_cancel", function (e) {
    var value = $(this).data('value');
    var td = $(this).parent("td");
    $(this).remove();
    $(this).closest('button').remove();
    td.html(value);
    console.log(value)
    td.addClass("editablerow_stock");


  });

  const csrftoken = getCookie("csrftoken");
  const sendToServer = (data, td) => {
    console.log(data);
    $.ajax({
      url: `/organism/updateStock/${data['name']}`,
      type: "POST",
      headers: { "X-CSRFToken": csrftoken },
      data: data,
    })
      .done((response) => {
        console.log(response);
      })
      .fail(() => {
        console.log("Error occured");
      });
  };
});
