$(document).ready(function () {
  console.log("loading data vitek");
  var a = "{{defaultvalues}}";
  console.log(a);
  $(".submit_data").click(function () {
    //selected objects
    console.log(a);

    var selected_data = [];

    $("input:checkbox[name=type]:checked").each(function () {
      selected_data.push($(this).val().toString());
    });
    var card_barcode = $("input[name=card_barcode]").val();
    // plot value
    var value_str = $("[data-name=data_process_value] option:selected")
      .val()
      .toString();
    console.log(value_str);
    // plot functijon
    var data_function_str = $("[data-name=data_function_type] option:selected")
      .val()
      .toString();

    // index value
    var index_value = [];
    $("#sortable3 li").each(function () {
      index_value.push($(this).text());
    });
    var index_values_str = index_value.toString();
    // column value
    var column_value = [];
    $("#sortable2 li").each(function () {
      column_value.push($(this).text());
    });
    var column_value_str = column_value.toString();
    // all data
    var data = {
      // data_map: data_map_str,

      selected_data: selected_data,
      values: value_str,
      columns: column_value_str,
      index: index_values_str,
      card_barcode: card_barcode,
      functions: data_function_str,
    };
    // console.log(data);
    // ajax send data to server, receive data from server
    sendToServer(data);
    // var table = result ? result : null;
    // console.log(table);
    // $("#pivotable").html = "";
    // $("#pivotable").html += table;
  });
});

const csrftoken = getCookie("csrftoken");
const sendToServer = (data) => {
  console.log("send to server");
  // console.log(data);

  $.ajax({
    url: "/vitekcard_list", //url,
    type: "POST",
    headers: { "X-CSRFToken": csrftoken },
    data: data,
  })
    .done((response) => {
      $("#pivotable").html("");
      if (response["msg"]) {
        saveData(response["table"], "pivottable.csv");
        $("#pivotable").append(response["msg"])
      } else {
        data = response["table"];
        // console.log(data);
        $("#pivotable").append(data);
      }
    })
    .fail((XMLHttpRequest, textStatus, errorThrown) => {
      console.log(XMLHttpRequest, textStatus, errorThrown);
    });
};
