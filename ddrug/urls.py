from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path
from rest_framework_simplejwt.views import TokenObtainPairView, TokenRefreshView


from ddrug.views import  (DrugListView, DrugCardView,detailDrug, DrugCreateView, DrugUpdateView,#createDrug, updateDrug, 
    VitekCard_ListView, VitekAST_ListView, VitekID_ListView, 
    smartsQuery, 
    ketcher_test,iframe_url, API_VITEK_ASTList, API_Drug_List, API_Drug_Detail, MIC_COADDListView, MIC_COADDCardView, 
    MIC_PubListView, MIC_PubCardView, MIC_PubListView, MIC_PubCardView, BreakpointListView)
from ddrug.upload_views import Import_VitekView, Import_DrugView


urlpatterns = [
    # API path
    path('api/token/', TokenObtainPairView.as_view(), name='token_obtain_pair'),
    path('api/token/refresh/', TokenRefreshView.as_view(), name='token_refresh'),
    path('api-vitek-ast/', API_VITEK_ASTList.as_view(), name="api_vitekast"),
    path('api-drug/', API_Drug_List.as_view(), name="api_drug"),
    path('api-drug/<str:pk>', API_Drug_Detail.as_view(), name="api_drug_detail"),

    # Normal path
    path('drug_card', DrugCardView.as_view(), name="drug_card"),
    path('drug_list', DrugListView.as_view(), name="drug_list"),
    path('drug/<str:pk>', detailDrug, name="drug_detail"),
    path('drug_detail_structure/<str:pk>', smartsQuery, name="smartsquery"),
    path('createDrug/', DrugCreateView.as_view(), name="drug_create"),
    path('updateDrug/<str:pk>', DrugUpdateView.as_view(), name="drug_update"),

    path('vitekcard_list', VitekCard_ListView.as_view(), name="vitekcard_list"),
    #path('vitekcard_detail/<str:pk>', detailVitekcard, name="vitekcard_detail"),
    path('vitekast_list', VitekAST_ListView.as_view(), name="vitekast_list"),
    path('vitekid_list', VitekID_ListView.as_view(), name="vitekid_list"),

    path('mic-coadd_list', MIC_COADDListView.as_view(), name="mic_coadd_list"),
    path('mic-coadd_card', MIC_COADDCardView.as_view(), name="mic_coadd_card"),
    path('mic-pub_list', MIC_PubListView.as_view(), name="mic_pub_list"),
    path('mic-pub_card', MIC_PubCardView.as_view(), name="mic_pub_card"),
    path('breakpoint', BreakpointListView.as_view(), name="breakpoint_list"),
     

    # path("import/<str:process_name>/", Importhandler_VITEK.as_view(), name="import-VITEK"),
    path("ketcher_test/", ketcher_test, name="ketcher_test"),

    path('import-vitek/', Import_VitekView.as_view(), name='import-vitek'),
    path('import-drug/', Import_DrugView.as_view(), name='import-drug'),

    
]