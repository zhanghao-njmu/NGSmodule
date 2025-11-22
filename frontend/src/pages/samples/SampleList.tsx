/**
 * Sample List Page - Complete CRUD for sample management
 */
import React, { useEffect, useState } from 'react'
import {
  Button,
  Space,
  Select,
  Upload,
  Tag,
  Modal,
  Form,
  Input,
  Popconfirm,
  Tooltip,
} from 'antd'
import {
  PlusOutlined,
  UploadOutlined,
  ExperimentOutlined,
  EditOutlined,
  DeleteOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useSampleStore } from '../../store/sampleStore'
import { useProjectStore } from '../../store/projectStore'
import { PageHeader, DataTable } from '../../components/common'
import { toast, notifications } from '../../utils/notification'
import type { Sample, SampleCreate, SampleUpdate } from '../../types/sample'
import dayjs from 'dayjs'

const { Option } = Select

export const SampleList: React.FC = () => {
  const { samples, loading, fetchSamples, createSample, updateSample, deleteSample, importFromCSV } =
    useSampleStore()
  const { projects, fetchProjects } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')
  const [isModalVisible, setIsModalVisible] = useState(false)
  const [editingSample, setEditingSample] = useState<Sample | null>(null)
  const [form] = Form.useForm()

  useEffect(() => {
    fetchProjects()
  }, [])

  useEffect(() => {
    if (selectedProject) {
      fetchSamples({ project_id: selectedProject })
    }
  }, [selectedProject])

  // 打开创建/编辑模态框
  const showModal = (sample?: Sample) => {
    if (sample) {
      setEditingSample(sample)
      form.setFieldsValue({
        sample_id: sample.sample_id,
        run_id: sample.run_id,
        group_name: sample.group_name,
        layout: sample.layout,
        batch_id: sample.batch_id,
      })
    } else {
      setEditingSample(null)
      form.resetFields()
    }
    setIsModalVisible(true)
  }

  // 关闭模态框
  const handleCancel = () => {
    setIsModalVisible(false)
    setEditingSample(null)
    form.resetFields()
  }

  // 提交表单（创建或更新）
  const handleSubmit = async () => {
    try {
      const values = await form.validateFields()

      if (editingSample) {
        // 更新样本
        const updateData: SampleUpdate = {
          sample_id: values.sample_id,
          run_id: values.run_id,
          group_name: values.group_name,
          layout: values.layout,
          batch_id: values.batch_id,
        }
        await updateSample(editingSample.id, updateData)
        toast.success('Sample updated successfully')
      } else {
        // 创建样本
        if (!selectedProject) {
          toast.warning('Please select a project first')
          return
        }
        const createData: SampleCreate = {
          project_id: selectedProject,
          sample_id: values.sample_id,
          run_id: values.run_id,
          group_name: values.group_name,
          layout: values.layout,
          batch_id: values.batch_id,
        }
        await createSample(createData)
        toast.success('Sample created successfully')
      }

      handleCancel()
      // 刷新列表
      if (selectedProject) {
        fetchSamples({ project_id: selectedProject })
      }
    } catch (error) {
      console.error('Failed to save sample:', error)
      toast.error('Failed to save sample')
    }
  }

  // 删除样本
  const handleDelete = async (id: string) => {
    try {
      await deleteSample(id)
      toast.success('Sample deleted successfully')
      // 刷新列表
      if (selectedProject) {
        fetchSamples({ project_id: selectedProject })
      }
    } catch (error) {
      console.error('Failed to delete sample:', error)
      toast.error('Failed to delete sample')
    }
  }

  // CSV 导入
  const handleCSVImport = async (file: File) => {
    if (!selectedProject) {
      toast.warning('Please select a project first')
      return false
    }

    const loadingToast = toast.loading('Importing samples...')
    try {
      await importFromCSV(selectedProject, file)
      loadingToast()
      notifications.uploadSuccess()
      // 刷新列表
      fetchSamples({ project_id: selectedProject })
    } catch (error) {
      loadingToast()
      notifications.uploadError()
    }
    return false // 阻止自动上传
  }

  const columns: ColumnsType<Sample> = [
    {
      title: 'Sample ID',
      dataIndex: 'sample_id',
      key: 'sample_id',
      render: (sample_id) => (
        <Space>
          <ExperimentOutlined style={{ color: 'var(--color-primary)' }} />
          <span style={{ fontWeight: 500 }}>{sample_id}</span>
        </Space>
      ),
    },
    {
      title: 'Run ID',
      dataIndex: 'run_id',
      key: 'run_id',
      width: 120,
      render: (run_id) => run_id || '-',
    },
    {
      title: 'Group',
      dataIndex: 'group_name',
      key: 'group_name',
      width: 120,
      render: (group) =>
        group ? (
          <Tag color={group.toLowerCase() === 'control' ? 'blue' : 'green'}>{group}</Tag>
        ) : (
          '-'
        ),
    },
    {
      title: 'Layout',
      dataIndex: 'layout',
      key: 'layout',
      width: 100,
      render: (layout) =>
        layout ? <Tag color={layout === 'PE' ? 'purple' : 'orange'}>{layout}</Tag> : '-',
    },
    {
      title: 'Batch',
      dataIndex: 'batch_id',
      key: 'batch_id',
      width: 120,
      render: (batch_id) => batch_id || '-',
    },
    {
      title: 'Files',
      dataIndex: 'file_count',
      key: 'file_count',
      width: 80,
      align: 'center',
      render: (count) => count || 0,
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).format('YYYY-MM-DD HH:mm'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 120,
      fixed: 'right',
      render: (_, record) => (
        <Space>
          <Tooltip title="Edit">
            <Button
              type="text"
              size="small"
              icon={<EditOutlined />}
              onClick={() => showModal(record)}
            />
          </Tooltip>
          <Popconfirm
            title="Delete Sample"
            description="Are you sure you want to delete this sample? All associated files will also be deleted."
            onConfirm={() => handleDelete(record.id)}
            okText="Yes, Delete"
            cancelText="Cancel"
            okButtonProps={{ danger: true }}
          >
            <Tooltip title="Delete">
              <Button type="text" size="small" danger icon={<DeleteOutlined />} />
            </Tooltip>
          </Popconfirm>
        </Space>
      ),
    },
  ]

  return (
    <div>
      <PageHeader
        left={
          <Select
            placeholder="Select Project"
            style={{ width: 300 }}
            value={selectedProject || undefined}
            onChange={setSelectedProject}
            loading={projects.length === 0}
          >
            {projects.map((p) => (
              <Option key={p.id} value={p.id}>
                {p.name}
              </Option>
            ))}
          </Select>
        }
        right={
          <>
            <Upload beforeUpload={handleCSVImport} accept=".csv" showUploadList={false}>
              <Button icon={<UploadOutlined />} disabled={!selectedProject}>
                Import CSV
              </Button>
            </Upload>
            <Button
              type="primary"
              icon={<PlusOutlined />}
              onClick={() => showModal()}
              disabled={!selectedProject}
            >
              New Sample
            </Button>
          </>
        }
      />

      <DataTable
        columns={columns}
        dataSource={samples}
        rowKey="id"
        loading={loading}
        pagination={{ pageSize: 20 }}
        emptyText="No Samples"
        emptyDescription="Select a project and add samples or import from CSV"
      />

      {/* 创建/编辑样本模态框 */}
      <Modal
        title={editingSample ? 'Edit Sample' : 'Create Sample'}
        open={isModalVisible}
        onOk={handleSubmit}
        onCancel={handleCancel}
        okText={editingSample ? 'Update' : 'Create'}
        width={600}
      >
        <Form form={form} layout="vertical">
          <Form.Item
            name="sample_id"
            label="Sample ID"
            rules={[{ required: true, message: 'Please enter sample ID' }]}
          >
            <Input placeholder="e.g., Sample001" />
          </Form.Item>

          <Form.Item name="run_id" label="Run ID">
            <Input placeholder="e.g., Run001 (optional)" />
          </Form.Item>

          <Form.Item name="group_name" label="Group Name">
            <Select placeholder="Select group (optional)" allowClear>
              <Option value="Control">Control</Option>
              <Option value="Treatment">Treatment</Option>
              <Option value="Case">Case</Option>
              <Option value="Healthy">Healthy</Option>
            </Select>
          </Form.Item>

          <Form.Item name="layout" label="Sequencing Layout">
            <Select placeholder="Select layout (optional)" allowClear>
              <Option value="PE">Paired-End (PE)</Option>
              <Option value="SE">Single-End (SE)</Option>
            </Select>
          </Form.Item>

          <Form.Item name="batch_id" label="Batch ID">
            <Input placeholder="e.g., Batch001 (optional)" />
          </Form.Item>
        </Form>
      </Modal>
    </div>
  )
}

export default SampleList
